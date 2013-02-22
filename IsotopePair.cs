using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.IO;

namespace Coon.NeuQuant
{
    class IsotopePair
    {
        public Pair parent;
        public int isotope;
        public double peak1MZ;
        public double peak1Intensity;
        public double peak2MZ;
        public double peak2Intensity;
        public double peak3MZ;
        public double peak3Intensity;
        public double peak4MZ;
        public double peak4Intensity;
        public double intensity;
        public bool missingChannel;

        // 2 isotopologue pairs
        public IsotopePair(Pair parent, int isotope, double lightMZ, double lightIntensity, double heavyMZ, double heavyIntensity)
        {
            this.parent = parent;
            this.isotope = isotope;
            this.peak1MZ = lightMZ;
            this.peak1Intensity = lightIntensity;
            this.peak2MZ = heavyMZ;
            this.peak2Intensity = heavyIntensity;
            if (lightMZ > 0 && heavyMZ > 0)
            {
                missingChannel = false;
                if (lightIntensity >= heavyIntensity)
                {
                    intensity = lightIntensity;
                }
                else
                {
                    intensity = heavyIntensity;
                }
            }
            else
            {
                missingChannel = true;
                if (lightMZ > 0)
                {
                    intensity = lightIntensity;
                }
                else
                {
                    intensity = heavyIntensity;
                }
            }

        }

        // 4 isotopologue pairs
        public IsotopePair(Pair parent, int isotope, double peak1MZ, double peak1Intensity, double peak2MZ, double peak2Intensity, double peak3MZ, double peak3Intensity, double peak4MZ, double peak4Intensity)
        {
            this.parent = parent;
            this.isotope = isotope;
            this.peak1MZ = peak1MZ;
            this.peak1Intensity = peak1Intensity;
            this.peak2MZ = peak2MZ;
            this.peak2Intensity = peak2Intensity;
            this.peak3MZ = peak3MZ;
            this.peak3Intensity = peak3Intensity;
            this.peak4MZ = peak4MZ;
            this.peak4Intensity = peak4Intensity;
            if (peak1MZ > 0 && peak2MZ > 0 && peak3MZ > 0 && peak4MZ > 0)
            {
                missingChannel = false;
            }
            else
            {
                missingChannel = true;
            }

            if (peak1Intensity > intensity)
            {
                intensity = peak1Intensity;
            }
            if (peak2Intensity > intensity)
            {
                intensity = peak2Intensity;
            }
            if (peak3Intensity > intensity)
            {
                intensity = peak3Intensity;
            }
            if (peak4Intensity > intensity)
            {
                intensity = peak4Intensity;
            }

        }
    }
}
