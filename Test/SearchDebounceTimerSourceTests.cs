using System.IO;
using System.Runtime.CompilerServices;
using NUnit.Framework;

namespace ProteaseGuru.Test
{
    /// <summary>
    /// Source-level regression guard for the shared search-debounce timer bug fixed on
    /// fix/protein-results-live-search-shared-timer.
    ///
    /// History: SearchModifications used to expose a single `public static DispatcherTimer Timer`.
    /// ProteinResultsWindow subscribed to it first; when IndividualProteinAnalyzerWindow became a
    /// second consumer (commit 34572d2, "Seek Maximum Coverage of One Protein by Protease Pairs and
    /// Triplets (#63)", 2026-04-06), its constructor called the same static SetUp(), which replaced
    /// the shared Timer object outright and re-subscribed Tick to itself -- so whichever window was
    /// constructed most recently silently "owned" every other window's debounced keystrokes. User
    /// report recorded 2026-09-21.
    ///
    /// Why source-level, not a WPF unit test: ProteaseGuruGui targets net10.0-windows with
    /// UseWPF=true, and DispatcherTimer requires the Windows desktop runtime to construct/tick. This
    /// (ProteaseGuru.Test) project is deliberately cross-platform (net10.0, no WPF reference) so it
    /// runs outside Windows too; adding a WPF project reference here would break that. These tests
    /// instead assert the fix's invariants directly against the source text -- exact, fast, and no
    /// UI thread required. The tradeoff: behavioral/dispatcher-level exercise of the fix (actually
    /// pumping the timer and confirming ticks route to the right window) still needs a person
    /// running the app on Windows; this Linux Gateway cannot execute WPF/DispatcherTimer code.
    /// </summary>
    public class SearchDebounceTimerSourceTests
    {
        private static string RepoRoot([CallerFilePath] string thisFilePath = "")
            // This file lives in <repoRoot>/Test/, so the repo root is one directory up.
            => Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFilePath)!, ".."));

        private static string ReadGuiSource(string relativePath)
        {
            var path = Path.Combine(RepoRoot(), "ProteaseGuruGui", relativePath);
            Assert.That(File.Exists(path), Is.True, $"Expected source file not found: {path}");
            return File.ReadAllText(path);
        }

        [Test]
        public void SearchModifications_ExposesNoStaticTimerOrSetUp()
        {
            var source = ReadGuiSource("SearchModifications.cs");

            Assert.That(source, Does.Not.Contain("static DispatcherTimer"),
                "SearchModifications must not go back to a static DispatcherTimer field; a static " +
                "timer is exactly the bug that let one window's constructor steal another window's " +
                "Tick subscription.");
            Assert.That(source, Does.Not.Contain("static void SetUp"),
                "SearchModifications must not reintroduce a static SetUp() factory method; timer " +
                "construction must stay on the instance constructor.");
        }

        [TestCase("ProteinResultsWindow.xaml.cs")]
        [TestCase("IndividualProteinAnalyzerWindow.xaml.cs")]
        public void ConsumingWindow_HoldsSearchDebounceAsAPrivateInstanceField(string fileName)
        {
            var source = ReadGuiSource(fileName);

            Assert.That(source, Does.Contain("private readonly SearchModifications _searchDebounce"),
                $"{fileName} must hold the debounce timer as its own private instance field, " +
                "constructed per window, so it can never be shared with another window.");
            Assert.That(source, Does.Not.Contain("SearchModifications.SetUp()"),
                $"{fileName} must not call a static SearchModifications.SetUp(); that call pattern " +
                "is what let a later window replace an earlier window's shared timer.");
            Assert.That(source, Does.Not.Contain("SearchModifications.Timer"),
                $"{fileName} must not reference a static SearchModifications.Timer; use the " +
                "instance's own _searchDebounce.Timer field instead.");
        }
    }
}
