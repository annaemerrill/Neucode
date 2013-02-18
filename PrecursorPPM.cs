using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon.Spectra;
using Coon;

namespace OMNE
{
    class PrecursorPPM : IComparable<PrecursorPPM>
    {
        public int Charge { get; set; }
        public Peptide Peptide { get; set; }
        public double EValue { get; set; }
        public double Ppm { get; set; }

        public PrecursorPPM(int charge, Peptide peptide, double eValue, double ppm)
        {
            Charge = charge;
            Peptide = peptide;
            EValue = eValue;
            Ppm = ppm;
        }

        public int CompareTo(PrecursorPPM other)
        {
            int comp = this.Ppm.CompareTo(other.Ppm);
            return comp;
        }

    }
}
