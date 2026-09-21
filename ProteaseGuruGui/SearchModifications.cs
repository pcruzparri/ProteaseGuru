using System;
using System.Windows.Threading;

namespace ProteaseGuru.Gui
{
    /// <summary>
    /// Debounces rapid search-box keystrokes behind a single DispatcherTimer tick.
    /// This type is intentionally instance-scoped: each owning window/control must
    /// construct and hold its own instance. A shared/static instance let a second
    /// window's constructor replace the first window's timer and steal its Tick
    /// subscription, so the first window's searches silently fired the second
    /// window's filter (or nothing at all) once both windows existed together
    /// (regression introduced in commit 34572d2, 2026-04-06, when
    /// IndividualProteinAnalyzerWindow became a second consumer of what was then
    /// a shared static timer).
    /// </summary>
    class SearchModifications
    {
        public DispatcherTimer Timer { get; }

        public SearchModifications()
        {
            Timer = new DispatcherTimer();
            Timer.Interval = TimeSpan.FromMilliseconds(300);
        }

        // starts timer to keep track of user keystrokes
        public void SetTimer()
        {
            // Reset the timer
            Timer.Stop();
            Timer.Start();
        }
    }
}
