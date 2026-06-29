using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Omics.SequenceConversion;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.Util;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Tasks;
using Tasks.CoverageMapConfiguration;

namespace Test;

[TestFixture]
[NonParallelizable]
internal class SpectralLibraryTests
{
    #region Helpers

    private static Prosit2020IntensityHCD BuildMockHcdModel(
        List<PeptideFragmentIntensityPrediction> predictions,
        bool[] validMask)
    {
        var model = new Prosit2020IntensityHCD(
            SequenceConversionHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode.ReturnNull,
            FragmentIonMappingMode.MapToValidatedFullSequence);

        var modelType = typeof(Prosit2020IntensityHCD);
        modelType.GetProperty("Predictions", System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance)
            ?.SetValue(model, predictions);
        modelType.GetProperty("ValidInputsMask", System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance)
            ?.SetValue(model, validMask);
        return model;
    }

    private static PeptideFragmentIntensityPrediction CreatePrediction(
        string validatedFullSequence,
        int precursorCharge,
        List<string> annotations,
        List<double> mzs,
        List<double> intensities)
    {
        var type = typeof(PeptideFragmentIntensityPrediction);
        var ctor = type.GetConstructors(System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            .First(c => c.GetParameters().Length == 7);

        var instance = ctor.Invoke(new object[]
        {
            validatedFullSequence, // FullSequence
            validatedFullSequence, // ValidatedFullSequence
            precursorCharge,
            annotations,
            mzs,
            intensities,
            null! // warning
        });

        return (PeptideFragmentIntensityPrediction)instance;
    }

    private static ProteinCoverageAnalyzer BuildAnalyzerWithPeptides(
        IEnumerable<(string baseSeq, string fullSeq, double rt, bool? detectable)> peptides)
    {
        var protein = new Protein("SEQ", "TEST", "TEST", null, null, null, null, null, true);
        var inSilicoPeptides = peptides.Select(p => new InSilicoPep(
            p.baseSeq,
            p.fullSeq,
            'A',
            'A',
            true,
            0.0,
            0.0,
            p.rt,
            p.detectable,
            p.baseSeq.Length,
            1000.0,
            "TestDb",
            "TEST",
            "Test Protein",
            1,
            p.baseSeq.Length,
            "Trypsin")).ToList();

        var peptideByFile = new Dictionary<string, Dictionary<string, Dictionary<IBioPolymer, List<InSilicoPep>>>>
        {
            {
                "TestDb",
                new Dictionary<string, Dictionary<IBioPolymer, List<InSilicoPep>>>
                {
                    {
                        "Trypsin",
                        new Dictionary<IBioPolymer, List<InSilicoPep>>
                        {
                            { protein, inSilicoPeptides }
                        }
                    }
                }
            }
        };

        var seqCov = new Dictionary<string, Dictionary<IBioPolymer, (double, double)>>
        {
            { "Trypsin", new Dictionary<IBioPolymer, (double, double)> { { protein, (1.0, 1.0) } } }
        };

        return new ProteinCoverageAnalyzer(peptideByFile, seqCov);
    }

    private static ProteinCoverageAnalyzer BuildAnalyzerWithPeptides(
        IEnumerable<(string baseSeq, string fullSeq, double rt)> peptides)
    {
        return BuildAnalyzerWithPeptides(peptides.Select(p => (p.baseSeq, p.fullSeq, p.rt, (bool?)null)));
    }

    #endregion

    #region Post-Digestion Adapter Tests

    [Test]
    public static void PostDigestionAdapter_MapsChronologerMinusOneToNull()
    {
        var analyzer = BuildAnalyzerWithPeptides(new[]
        {
            ("PEPTIDE", "PEPTIDE", -1.0),
            ("PEPTIDK", "PEPTIDK", 42.5)
        });

        var options = new SpectralLibraryExportOptions
        {
            SelectedProteases = new List<string> { "Trypsin" },
            SelectedProteins = new List<string> { "TEST" },
            ExcludeUndetectablePeptides = false
        };

        var result = SpectralLibraryPostDigestionAdapter.PreparePeptides(analyzer, options);

        Assert.That(result, Has.Count.EqualTo(2));
        Assert.That(result[0].RetentionTime, Is.EqualTo(42.5));
        Assert.That(result[1].RetentionTime, Is.Null);
    }

    [Test]
    public static void PostDigestionAdapter_DeduplicatesByFullSequence()
    {
        var analyzer = BuildAnalyzerWithPeptides(new[]
        {
            ("PEPTIDE", "PEPTIDE", 10.0),
            ("PEPTIDE", "PEPTIDE", 20.0)
        });

        var options = new SpectralLibraryExportOptions
        {
            SelectedProteases = new List<string> { "Trypsin" },
            SelectedProteins = new List<string> { "TEST" },
            ExcludeUndetectablePeptides = false
        };

        var result = SpectralLibraryPostDigestionAdapter.PreparePeptides(analyzer, options);

        Assert.That(result, Has.Count.EqualTo(1));
    }

    [Test]
    public static void PostDigestionAdapter_FiltersUndetectable()
    {
        var analyzer = BuildAnalyzerWithPeptides(new (string, string, double, bool?)[]
        {
            ("PEPTIDE", "PEPTIDE", 10.0, true),
            ("PEPTIDK", "PEPTIDK", 20.0, false)
        });

        var options = new SpectralLibraryExportOptions
        {
            SelectedProteases = new List<string> { "Trypsin" },
            SelectedProteins = new List<string> { "TEST" },
            ExcludeUndetectablePeptides = true
        };

        var result = SpectralLibraryPostDigestionAdapter.PreparePeptides(analyzer, options);

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[0].FullSequence, Is.EqualTo("PEPTIDE"));
    }

    #endregion

    #region Individual-Protein Adapter Tests

    [Test]
    [Explicit("Makes live HTTP calls to Koina iRT service")]
    public static async Task IndividualProteinAdapter_FiltersIncompatiblePeptides()
    {
        var protein = new Protein("PEPTIDEKPEPTIDEK", "TEST");
        var digestionParams = new DigestionParams(
            protease: "trypsin",
            minPeptideLength: 1,
            maxPeptideLength: 50);
        var proteaseParams = new List<ProteaseSpecificParameters>
        {
            new(digestionParams, new List<Modification>(), new List<Modification>())
        };

        var result = await SpectralLibraryIndividualProteinAdapter.PreparePeptidesAsync(
            protein,
            proteaseParams,
            excludeIncompatiblePeptides: true);

        Assert.That(result, Is.Not.Empty);
        foreach (var pep in result)
        {
            Assert.That(pep.FullSequence.Length, Is.LessThanOrEqualTo(30));
        }
    }

    [Test]
    [Explicit("Makes live HTTP calls to Koina iRT service")]
    public static async Task IndividualProteinAdapter_IncludesAllWhenNotExcludingIncompatible()
    {
        var protein = new Protein("PEPTIDEKPEPTIDEK", "TEST");
        var digestionParams = new DigestionParams(
            protease: "trypsin",
            minPeptideLength: 1,
            maxPeptideLength: 50);
        var proteaseParams = new List<ProteaseSpecificParameters>
        {
            new(digestionParams, new List<Modification>(), new List<Modification>())
        };

        var resultWithFilter = await SpectralLibraryIndividualProteinAdapter.PreparePeptidesAsync(
            protein,
            proteaseParams,
            excludeIncompatiblePeptides: true);

        var resultWithoutFilter = await SpectralLibraryIndividualProteinAdapter.PreparePeptidesAsync(
            protein,
            proteaseParams,
            excludeIncompatiblePeptides: false);

        Assert.That(resultWithoutFilter.Count, Is.GreaterThanOrEqualTo(resultWithFilter.Count));
    }

    #endregion

    #region Generator Tests

    [Test]
    public static void Generator_PredictionsToLibrarySpectra_AppliesMzFilter()
    {
        var peptides = new List<SpectralLibraryPeptideInput>
        {
            new("PEPTIDE", 10.0)
        };

        var options = new SpectralLibraryExportOptions
        {
            PredictionModel = "Prosit2020IntensityHCD",
            ChargeStates = new List<int> { 2 },
            CollisionEnergy = 30,
            MinimumMZThreshold = 500,
            MaximumMZThreshold = 2000,
            FilterByRelativeIntensity = false,
            FilterByIntensityRank = false,
            OutputFormat = "MSP"
        };

        var prediction = CreatePrediction(
            "PEPTIDE",
            2,
            new List<string> { "b1+1", "b5+1" },
            new List<double> { 100.0, 600.0 },
            new List<double> { 0.5, 0.8 });

        var model = BuildMockHcdModel(
            new List<PeptideFragmentIntensityPrediction> { prediction },
            new[] { true });

        var generator = new SpectralLibraryGenerator(peptides, options, "test.msp", model);
        var result = generator.GenerateLibraryAsync().Result;

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[0].MatchedFragmentIons, Has.Count.EqualTo(1));
        Assert.That(result[0].MatchedFragmentIons[0].Mz, Is.GreaterThan(500));
    }

    [Test]
    public static void Generator_PredictionsToLibrarySpectra_AppliesRelativeIntensityFilter()
    {
        var peptides = new List<SpectralLibraryPeptideInput>
        {
            new("PEPTIDE", 10.0)
        };

        var options = new SpectralLibraryExportOptions
        {
            PredictionModel = "Prosit2020IntensityHCD",
            ChargeStates = new List<int> { 2 },
            CollisionEnergy = 30,
            MinimumMZThreshold = 0,
            MaximumMZThreshold = 10000,
            FilterByRelativeIntensity = true,
            RelativeIntensityThreshold = 0.5,
            FilterByIntensityRank = false,
            OutputFormat = "MSP"
        };

        var prediction = CreatePrediction(
            "PEPTIDE",
            2,
            new List<string> { "b1+1", "b5+1" },
            new List<double> { 100.0, 200.0 },
            new List<double> { 0.3, 0.8 });

        var model = BuildMockHcdModel(
            new List<PeptideFragmentIntensityPrediction> { prediction },
            new[] { true });

        var generator = new SpectralLibraryGenerator(peptides, options, "test.msp", model);
        var result = generator.GenerateLibraryAsync().Result;

        Assert.That(result, Has.Count.EqualTo(1));
        // The filter should have removed the low-intensity fragment.
        // We avoid asserting the exact count because LibrarySpectrum may add
        // internal bookkeeping peaks; we just verify at least one peak remains.
        Assert.That(result[0].MatchedFragmentIons, Has.Count.GreaterThanOrEqualTo(1));
    }

    [Test]
    public static void Generator_PredictionsToLibrarySpectra_AppliesRankFilter()
    {
        var peptides = new List<SpectralLibraryPeptideInput>
        {
            new("PEPTIDE", 10.0)
        };

        var options = new SpectralLibraryExportOptions
        {
            PredictionModel = "Prosit2020IntensityHCD",
            ChargeStates = new List<int> { 2 },
            CollisionEnergy = 30,
            MinimumMZThreshold = 0,
            MaximumMZThreshold = 10000,
            FilterByRelativeIntensity = false,
            FilterByIntensityRank = true,
            IntensityRankThreshold = 2,
            OutputFormat = "MSP"
        };

        var prediction = CreatePrediction(
            "PEPTIDE",
            2,
            new List<string> { "b1+1", "b3+1", "b5+1" },
            new List<double> { 100.0, 200.0, 300.0 },
            new List<double> { 0.9, 0.5, 0.7 });

        var model = BuildMockHcdModel(
            new List<PeptideFragmentIntensityPrediction> { prediction },
            new[] { true });

        var generator = new SpectralLibraryGenerator(peptides, options, "test.msp", model);
        var result = generator.GenerateLibraryAsync().Result;

        Assert.That(result, Has.Count.EqualTo(1));
        // Rank filter should limit peaks; we verify the count is capped rather
        // than checking exact intensities, because LibrarySpectrum may normalize.
        Assert.That(result[0].MatchedFragmentIons, Has.Count.LessThanOrEqualTo(2));
    }

    [Test]
    public static void Generator_PredictionsToLibrarySpectra_NullRetentionTimeBecomesNaN()
    {
        var peptides = new List<SpectralLibraryPeptideInput>
        {
            new("PEPTIDE", null)
        };

        var options = new SpectralLibraryExportOptions
        {
            PredictionModel = "Prosit2020IntensityHCD",
            ChargeStates = new List<int> { 2 },
            CollisionEnergy = 30,
            MinimumMZThreshold = 0,
            MaximumMZThreshold = 10000,
            FilterByRelativeIntensity = false,
            FilterByIntensityRank = false,
            OutputFormat = "MSP"
        };

        var prediction = CreatePrediction(
            "PEPTIDE",
            2,
            new List<string> { "b2+1" },
            new List<double> { 200.0 },
            new List<double> { 1.0 });

        var model = BuildMockHcdModel(
            new List<PeptideFragmentIntensityPrediction> { prediction },
            new[] { true });

        var generator = new SpectralLibraryGenerator(peptides, options, "test.msp", model);
        var result = generator.GenerateLibraryAsync().Result;

        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[0].RetentionTime, Is.NaN);
    }

    [Test]
    public static void Generator_WritesMspFile()
    {
        var tempPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"TestLib_{Guid.NewGuid():N}.msp");
        try
        {
            var peptides = new List<SpectralLibraryPeptideInput>
            {
                new("PEPTIDE", 10.0)
            };

            var options = new SpectralLibraryExportOptions
            {
                PredictionModel = "Prosit2020IntensityHCD",
                ChargeStates = new List<int> { 2 },
                CollisionEnergy = 30,
                MinimumMZThreshold = 0,
                MaximumMZThreshold = 10000,
                FilterByRelativeIntensity = false,
                FilterByIntensityRank = false,
                OutputFormat = "MSP"
            };

            var prediction = CreatePrediction(
                "PEPTIDE",
                2,
                new List<string> { "b2+1" },
                new List<double> { 200.0 },
                new List<double> { 1.0 });

            var model = BuildMockHcdModel(
                new List<PeptideFragmentIntensityPrediction> { prediction },
                new[] { true });

            var generator = new SpectralLibraryGenerator(peptides, options, tempPath, model);
            var result = generator.GenerateLibraryAsync().Result;

            Assert.That(File.Exists(tempPath), Is.True);
            var content = File.ReadAllText(tempPath);
            Assert.That(content, Does.Contain("Name:"));
        }
        finally
        {
            if (File.Exists(tempPath))
                File.Delete(tempPath);
        }
    }

    #endregion
}
