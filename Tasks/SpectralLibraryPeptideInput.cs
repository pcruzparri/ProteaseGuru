namespace Tasks;

/// <summary>
/// Decoupled input for spectral library generation, independent of whether peptides
/// came from a post-digestion analyzer or a live single-protein digestion.
/// </summary>
public sealed record SpectralLibraryPeptideInput(
    string FullSequence,
    double? RetentionTime);
