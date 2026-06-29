using Tasks.CoverageMapConfiguration;

namespace Tasks;

/// <summary>
/// Adapter that prepares peptide inputs for spectral-library generation from
/// post-digestion results (InSilicoPep objects produced by DigestionTask).
/// </summary>
public static class SpectralLibraryPostDigestionAdapter
{
    /// <summary>
    /// Gathers, filters, deduplicates, and maps peptides from the analyzer into the
    /// unified DTO used by <see cref="SpectralLibraryGenerator"/>.
    /// </summary>
    public static List<SpectralLibraryPeptideInput> PreparePeptides(
        ProteinCoverageAnalyzer analyzer,
        SpectralLibraryExportOptions options)
    {
        var peptides = new List<InSilicoPep>();

        var selectedProteins = analyzer.ProteinCoverageResults.Keys
            .Where(p => options.SelectedProteins.Contains(p.Accession))
            .ToList();

        foreach (var protein in selectedProteins)
        {
            foreach (var proteaseName in options.SelectedProteases)
            {
                var proteinPeptides = analyzer.GetPeptidesForProteinAndProtease(protein, proteaseName);
                peptides.AddRange(proteinPeptides);
            }
        }

        if (options.ExcludeUndetectablePeptides)
        {
            peptides = peptides.Where(p => p.PflyDetectability == true).ToList();
        }

        // Deduplicate by FullSequence. The spectral library is indexed by
        // sequence/charge, so modification variants with the same FullSequence
        // would collide anyway. Order deterministically so the RT value is stable.
        return peptides
            .OrderBy(p => p.ChronologerRetentionTime == -1 ? 1 : 0)
            .ThenBy(p => p.ChronologerRetentionTime)
            .DistinctBy(p => p.FullSequence)
            .Select(p => new SpectralLibraryPeptideInput(
                p.FullSequence,
                p.ChronologerRetentionTime == -1 ? null : p.ChronologerRetentionTime))
            .ToList();
    }
}
