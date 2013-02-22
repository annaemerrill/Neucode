using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL.Proteomics;
using CSMSL;

namespace Coon.NeuQuant
{
    class PrecursorPPM : IComparable<PrecursorPPM>
    {
        public int Charge { get; set; }
        public string Peptide { get; set; }
        public double EValue { get; set; }
        public double Ppm { get; set; }

        public PrecursorPPM(int charge, string peptide, double eValue, double ppm)
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
