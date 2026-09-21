using System.IO;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;
using NUnit.Framework;

namespace ProteaseGuru.Test
{
    /// <summary>
    /// Source-level guard that the search debounce timer stays instance-scoped per window.
    /// ProteaseGuruGui is net10.0-windows/WPF and DispatcherTimer needs the Windows desktop
    /// runtime, so it can't be exercised from this cross-platform test project; these tests
    /// assert the invariant against the source text instead. Behavioral verification (both
    /// windows' searches filtering independently) is manual on Windows.
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

            Assert.That(Regex.IsMatch(source,
                    @"\bstatic\b[^;{}\r\n]*\bDispatcherTimer\b|\bDispatcherTimer\b[^;{}\r\n]*\bstatic\b"),
                Is.False,
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

            Assert.That(Regex.IsMatch(source,
                    @"\bprivate\s+readonly\s+SearchModifications\s+_searchDebounce\s*=\s*new\s+SearchModifications\s*\(\s*\)\s*;"),
                Is.True,
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
