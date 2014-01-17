using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL.Spectral;
using CSMSL.Chemistry;
using CSMSL.IO.Thermo;
using CSMSL;

namespace Coon.NeuQuant
{
    public class Spacing
    {
        public string Sequence;
        public double RetentionTime;
        public int Charge;
        public int NumLabels;
        public int Isotope;
        public double LightInt;
        public double HeavyInt;
        public double LightMZ;
        public double HeavyMZ;
        public double TheoSpacingDa;
        public double TheoSpacingTh;
        public double SpacingDa;
        public double SpacingTh;
        public int ScanNumber;
        public double Log10SummedIntensity;
        public double SpacingMTh;

        public Spacing(string sequence, double retentionTime, int scanNumber, int charge, int numLabels, int isotope, PeptideID peptide, ILabeledPeak light, ILabeledPeak heavy = null)
        {
            Sequence = sequence;
            RetentionTime = retentionTime;
            Charge = charge;
            NumLabels = numLabels;
            Isotope = isotope;
            TheoSpacingDa = peptide.spacingMassRange[1,0].Mean;
            TheoSpacingTh = TheoSpacingDa / Charge;
            ScanNumber = scanNumber;

            if (light != null)
            {
                LightInt = light.Y;
                LightMZ = light.X;
            }

            if (heavy != null)
            {
                HeavyInt = heavy.Y;
                HeavyMZ = heavy.X;

                SpacingDa = Mass.MassFromMz(HeavyMZ, Charge) - Mass.MassFromMz(LightMZ, Charge);
                SpacingTh = HeavyMZ - LightMZ;
                SpacingMTh = SpacingTh * 1000.0;
                Log10SummedIntensity = Math.Log10(LightInt + HeavyInt);
            }
            else
            {
                HeavyInt = 0;
                HeavyMZ = 0;

                SpacingDa = 0.0;
                SpacingTh = 0.0;
            }
        }
    }
}
