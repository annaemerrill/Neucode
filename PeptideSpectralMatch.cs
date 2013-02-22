using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL;
using CSMSL.Proteomics;

namespace Coon.NeuQuant
{
    class PeptideSpectralMatch: IComparable<PeptideSpectralMatch>, IComparable<double>, IComparable<int>
    {
        public int ScanNumber;
        public int Charge;
        public double EValue;

        public PeptideSpectralMatch(int scanNumber, int charge, double eValue)
        {
            ScanNumber = scanNumber;
            Charge = charge;
            EValue = eValue;
        }

        public int CompareTo(PeptideSpectralMatch other)
        {
            if (other == null) return 1;

            return EValue.CompareTo(other.EValue);
        }

        public int CompareTo(double other)
        {
            if (other == double.NaN) return 1;

            return EValue.CompareTo(other);
        }

        public int CompareTo(int other)
        {
            if (other == 0) return 1;

            return ScanNumber.CompareTo(other);
        }
    }
}
