using System.Text.RegularExpressions;
using Omics;
using Omics.Modifications;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Proteomics.ProteolyticDigestion;

namespace Tasks;

/// <summary>
/// Adapter that prepares peptide inputs for spectral-library generation from a
/// single protein by performing a live digestion, filtering, and iRT prediction.
/// </summary>
public static class SpectralLibraryIndividualProteinAdapter
{
    private static readonly Regex ModPattern = new(@"\[[^\]]+\]", RegexOptions.Compiled);

    /// <summary>
    /// Digests the protein with the supplied protease parameters, filters to
    /// Prosit-compatible peptides, predicts iRT values, and maps to the unified DTO.
    /// </summary>
    public static async Task<List<SpectralLibraryPeptideInput>> PreparePeptidesAsync(
        IBioPolymer protein,
        IEnumerable<ProteaseSpecificParameters> proteaseParams,
        bool excludeIncompatiblePeptides,
        IProgress<string>? progress = null,
        CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();

        var uniqueMzLibSequences = DigestToUniqueMzLibSequences(
            protein, proteaseParams, excludeIncompatiblePeptides);

        if (uniqueMzLibSequences.Count == 0)
        {
            throw new InvalidOperationException(
                "No peptides passed the Prosit filter. " +
                "Try adjusting missed cleavages or peptide length limits.");
        }

        progress?.Report($"Digest complete: {uniqueMzLibSequences.Count} unique Prosit-compatible peptides.");

        // Predict iRT via Prosit2019iRT (Koina). Wrapped in Task.Run because
        // Predict() makes batched HTTP requests synchronously.
        progress?.Report("Requesting iRT values from Prosit_2019_irt via Koina...");
        cancellationToken.ThrowIfCancellationRequested();

        var rtModel = new Prosit2019iRT();
        var rtInputs = uniqueMzLibSequences
            .Select(s => new RetentionTimePredictionInput(s))
            .ToList();

        var rtPredictions = await Task.Run(() => rtModel.Predict(rtInputs), cancellationToken);

        int rtValidCount = rtPredictions.Count(p => p.PredictedRetentionTime.HasValue);
        if (rtValidCount == 0)
        {
            string firstReason = rtPredictions
                .Select(p => p.Warning?.Message)
                .FirstOrDefault(m => !string.IsNullOrEmpty(m)) ?? "no details";
            throw new InvalidOperationException(
                $"All peptides were rejected by Prosit_2019_irt. Details: {firstReason}");
        }

        var mzLibToIrt = new Dictionary<string, double?>();
        foreach (var pred in rtPredictions)
        {
            if (pred.FullSequence != null)
                mzLibToIrt[pred.FullSequence] = pred.PredictedRetentionTime;
        }

        progress?.Report($"iRT predictions received for {rtValidCount} peptides.");

        return uniqueMzLibSequences
            .Select(seq => new SpectralLibraryPeptideInput(
                seq,
                mzLibToIrt.TryGetValue(seq, out var rt) ? rt : null))
            .ToList();
    }

    /// <summary>
    /// Digests the protein with every supplied protease, deduplicates by FullSequence,
    /// and optionally pre-filters to peptides that Prosit will accept.
    /// </summary>
    private static List<string> DigestToUniqueMzLibSequences(
        IBioPolymer protein,
        IEnumerable<ProteaseSpecificParameters> proteaseParams,
        bool excludeIncompatiblePeptides)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal);
        var result = new List<string>();

        foreach (var pp in proteaseParams)
        {
            var peptides = protein.Digest(pp.DigestionParams, pp.FixedMods, pp.VariableMods);

            foreach (var pep in peptides)
            {
                if (!seen.Add(pep.FullSequence))
                    continue;

                if (excludeIncompatiblePeptides)
                {
                    if (!IsPrositCompatible(pep.FullSequence))
                        continue;
                }

                result.Add(pep.FullSequence);
            }
        }

        return result;
    }

    /// <summary>
    /// Returns true if the mzLib-format full sequence is compatible with
    /// Prosit2019iRT and Prosit2020IntensityHCD.
    /// </summary>
    private static bool IsPrositCompatible(string fullSequence)
    {
        string baseSeq = ModPattern.Replace(fullSequence, string.Empty);

        if (baseSeq.Length is < 1 or > 30)
            return false;

        if (!baseSeq.All(c => "ACDEFGHIKLMNPQRSTVWY".Contains(c)))
            return false;

        bool allModsOk = ModPattern
            .Matches(fullSequence)
            .All(m =>
                m.Value == "[Common Fixed:Carbamidomethyl on C]" ||
                m.Value == "[Common Variable:Oxidation on M]");

        return allModsOk;
    }
}
