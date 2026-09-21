using System;
using System.Windows.Threading;

namespace ProteaseGuru.Gui
{
    /// <summary>Per-window debounce timer for a search box; each window owns its own instance.</summary>
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
            Timer.Stop();
            Timer.Start();
        }
    }
}
