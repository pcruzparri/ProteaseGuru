using System.Diagnostics;
using Chemistry;
using Omics.Fragmentation;
using Omics.SequenceConversion;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;

namespace Tasks;

/// <summary>
/// Configuration options for spectral library generation.
/// </summary>
public class SpectralLibraryExportOptions
{
    // Peptide source filtering options
    public List<string> SelectedProteases { get; set; }
    public List<string> SelectedProteins { get; set; }

    // Prediction model options
    public string PredictionModel { get; set; }
    public List<int> ChargeStates { get; set; }
    public int CollisionEnergy { get; set; }

    // Peptide filtering options
    public bool ExcludeIncompatiblePeptides { get; set; }
    public bool ExcludeUndetectablePeptides { get; set; }

    // Fragment ion filtering options
    public double MinimumMZThreshold { get; set; }
    public double MaximumMZThreshold { get; set; }
    public bool FilterByRelativeIntensity { get; set; }
    public double RelativeIntensityThreshold { get; set; }
    public bool FilterByIntensityRank { get; set; }
    public int IntensityRankThreshold { get; set; }

    // Output options
    public string OutputFormat { get; set; }
}

/// <summary>
/// Generates a spectral library from peptide inputs by predicting fragment
/// intensities and writing the result to disk.
/// </summary>
public class SpectralLibraryGenerator
{
    private readonly List<SpectralLibraryPeptideInput> _peptides;
    private readonly SpectralLibraryExportOptions _options;
    private readonly string _outputPath;
    private readonly FragmentIntensityModel? _injectedModel;

    /// <summary>
    /// Production constructor. The model is built from <see cref="SpectralLibraryExportOptions.PredictionModel"/>.
    /// </summary>
    public SpectralLibraryGenerator(
        IEnumerable<SpectralLibraryPeptideInput> peptides,
        SpectralLibraryExportOptions options,
        string outputPath)
    {
        _peptides = peptides?.ToList() ?? throw new ArgumentNullException(nameof(peptides));
        _options = options ?? throw new ArgumentNullException(nameof(options));
        _outputPath = outputPath ?? throw new ArgumentNullException(nameof(outputPath));

        if (_options.ChargeStates == null || _options.ChargeStates.Count == 0)
            throw new ArgumentException("At least one charge state must be specified.", nameof(options.ChargeStates));
    }

    /// <summary>
    /// Test constructor that injects a pre-built <see cref="FragmentIntensityModel"/>
    /// so unit tests can avoid network calls.
    /// </summary>
    public SpectralLibraryGenerator(
        IEnumerable<SpectralLibraryPeptideInput> peptides,
        SpectralLibraryExportOptions options,
        string outputPath,
        FragmentIntensityModel model)
    {
        _peptides = peptides?.ToList() ?? throw new ArgumentNullException(nameof(peptides));
        _options = options ?? throw new ArgumentNullException(nameof(options));
        _outputPath = outputPath ?? throw new ArgumentNullException(nameof(outputPath));
        _injectedModel = model ?? throw new ArgumentNullException(nameof(model));
    }

    /// <summary>
    /// Generates the spectral library asynchronously. Network-bound prediction
    /// steps run on thread-pool threads.
    /// </summary>
    public async Task<List<LibrarySpectrum>> GenerateLibraryAsync(
        IProgress<string>? progress = null,
        CancellationToken cancellationToken = default)
    {
        progress?.Report("Building fragment-intensity model...");
        cancellationToken.ThrowIfCancellationRequested();

        FragmentIntensityModel model = _injectedModel ?? BuildModel();

        var inputs = new List<FragmentIntensityPredictionInput>();
        var rts = new List<double?>();

        foreach (int pc in _options.ChargeStates)
        {
            inputs.AddRange(_peptides.Select(p => new FragmentIntensityPredictionInput(
                FullSequence: p.FullSequence,
                PrecursorCharge: pc,
                CollisionEnergy: _options.CollisionEnergy,
                InstrumentType: null,
                FragmentationType: null)));
            rts.AddRange(_peptides.Select(p => p.RetentionTime));
        }

        if (_injectedModel == null)
        {
            progress?.Report(
                $"Requesting HCD intensities from {_options.PredictionModel} " +
                $"({inputs.Count} spectra, NCE {_options.CollisionEnergy})...");
            cancellationToken.ThrowIfCancellationRequested();

            await Task.Run(() => model.Predict(inputs), cancellationToken);
        }

        progress?.Report("Processing predictions and applying filters...");
        cancellationToken.ThrowIfCancellationRequested();

        var library = await Task.Run(() => PredictionsToLibrarySpectra(model, rts), cancellationToken);

        cancellationToken.ThrowIfCancellationRequested();
        progress?.Report($"Writing {library.Count} spectra to {_outputPath}...");
        WriteLibrary(library);

        progress?.Report($"Done. {library.Count} spectra written.");
        return library;
    }

    private FragmentIntensityModel BuildModel()
    {
        switch (_options.PredictionModel)
        {
            case "Prosit2020IntensityHCD":
                return new Prosit2020IntensityHCD(
                    modHandlingMode: _options.ExcludeIncompatiblePeptides
                        ? SequenceConversionHandlingMode.ReturnNull
                        : SequenceConversionHandlingMode.RemoveIncompatibleElements,
                    parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull,
                    fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence);
            default:
                throw new NotSupportedException(
                    $"Prediction model {_options.PredictionModel} is not supported.");
        }
    }

    /// <summary>
    /// Mirrors mzLib's <c>FragmentIntensityModel.GenerateLibrarySpectraFromPredictions</c>,
    /// but adds m/z range, relative-intensity, and top-N rank filters.
    /// </summary>
    private List<LibrarySpectrum> PredictionsToLibrarySpectra(
        FragmentIntensityModel model,
        List<double?> retentionTimes)
    {
        Debug.Assert(model.Predictions.Count == model.ValidInputsMask.Length,
            "Predictions are expected to be realigned to the full input length.");

        var validPredictions = model.ValidInputsMask
            .Select((isValid, index) => (isValid, index))
            .Where(x => x.isValid)
            .Select(x => (Prediction: model.Predictions[x.index], RetentionTime: retentionTimes[x.index]))
            .ToList();

        var predictedSpectra = new List<LibrarySpectrum>();

        foreach (var (prediction, retentionTime) in validPredictions)
        {
            var peptide = new PeptideWithSetModifications(prediction.ValidatedFullSequence);
            var fragmentIons = new List<MatchedFragmentIon>();

            var theoreticalProducts = new List<Product>();
            peptide.Fragment(MassSpectrometry.DissociationType.HCD, FragmentationTerminus.Both, theoreticalProducts);

            var tpLookup = theoreticalProducts.ToDictionary(tp => tp.Annotation);
            var predictionAnnotationIntensityLookup = new Dictionary<string, double>();
            double maxFragmentIntensity = prediction.FragmentIntensities.DefaultIfEmpty(0).Max();

            for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
            {
                if (prediction.FragmentAnnotations[i] == null ||
                    !prediction.FragmentAnnotations[i].Contains("+") ||
                    prediction.FragmentMZs[i] < _options.MinimumMZThreshold ||
                    prediction.FragmentMZs[i] > _options.MaximumMZThreshold ||
                    (_options.FilterByRelativeIntensity &&
                     prediction.FragmentIntensities[i] < maxFragmentIntensity * _options.RelativeIntensityThreshold))
                {
                    continue;
                }
                predictionAnnotationIntensityLookup[prediction.FragmentAnnotations[i]] = prediction.FragmentIntensities[i];
            }

            foreach (var pa in predictionAnnotationIntensityLookup.Keys)
            {
                var productTypeAndCharge = pa.Split("+");
                if (productTypeAndCharge.Length < 2 ||
                    !tpLookup.TryGetValue(productTypeAndCharge[0], out var tp) ||
                    !int.TryParse(productTypeAndCharge[1], out int charge))
                {
                    continue;
                }

                var fragmentIon = new MatchedFragmentIon(
                    neutralTheoreticalProduct: tp,
                    experMz: tp.ToMz(charge),
                    experIntensity: predictionAnnotationIntensityLookup[pa],
                    charge: charge);

                fragmentIons.Add(fragmentIon);
            }

            if (_options.FilterByIntensityRank && _options.IntensityRankThreshold != -1)
            {
                fragmentIons = fragmentIons
                    .OrderByDescending(fi => fi.Intensity)
                    .Take(_options.IntensityRankThreshold)
                    .ToList();
            }

            var spectrum = new LibrarySpectrum(
                sequence: peptide.FullSequence,
                precursorMz: peptide.ToMz(prediction.PrecursorCharge),
                chargeState: prediction.PrecursorCharge,
                peaks: fragmentIons,
                rt: retentionTime ?? double.NaN);

            predictedSpectra.Add(spectrum);
        }

        // LibrarySpectrum.Name is "Sequence/ChargeState", so this only collapses genuine duplicates.
        return predictedSpectra.DistinctBy(p => p.Name).ToList();
    }

    private void WriteLibrary(List<LibrarySpectrum> spectra)
    {
        switch (_options.OutputFormat)
        {
            case "MSP":
                WriteMSP(spectra);
                break;
            default:
                throw new NotSupportedException(
                    $"Output format {_options.OutputFormat} is not supported.");
        }
    }

    private void WriteMSP(List<LibrarySpectrum> spectra)
    {
        var spectralLibrary = new SpectralLibrary();
        spectralLibrary.Results = spectra;
        spectralLibrary.WriteResults(_outputPath);
    }
}
