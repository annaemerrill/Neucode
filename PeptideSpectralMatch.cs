using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL;
using CSMSL.Proteomics;
using CSMSL.IO;

namespace Coon.NeuQuant
{
    public class PeptideSpectralMatch: IComparable<PeptideSpectralMatch>, IComparable<double>, IComparable<int>
    {
        public int ScanNumber;
        public int MS1ScanNumber;
        public int Charge;
        public double EValue;
        public double Resolvability;
        public MSDataFile RawFile;
        public string FilenameID;

        public PeptideSpectralMatch(int scanNumber, MSDataFile rawFile, int charge, double eValue, string fileNameID, double resolvability = double.NaN)
        {
            ScanNumber = scanNumber;
            Charge = charge;
            EValue = eValue;
            Resolvability = resolvability;
            RawFile = rawFile;
            FilenameID = fileNameID;

            bool stop = false;
            int currentScanNumber = ScanNumber;
            rawFile.Open();
            while (!stop && currentScanNumber >= rawFile.FirstSpectrumNumber)
            {
                if (rawFile[currentScanNumber].MsnOrder == 1) stop = true;
                else currentScanNumber--;
            }
            MS1ScanNumber = currentScanNumber;
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
