using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Coon.NeuQuant
{
    class CoalescenceCheck : IComparable<CoalescenceCheck>
    {
        public double intensity;
        public double missingChannelFrequency;

        public CoalescenceCheck(double Intensity, double MissingChannelFrequency)
        {
            intensity = Intensity;
            missingChannelFrequency = MissingChannelFrequency;
        }

        public int CompareTo(CoalescenceCheck other)
        {
            int comp = this.intensity.CompareTo(other.intensity);
            return comp;
        }
    }
}
