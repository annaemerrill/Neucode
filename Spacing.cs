using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL.Spectral;
using CSMSL.IO.Thermo;
using CSMSL;

namespace Coon.NeuQuant
{
    class Spacing : IComparable<Spacing>
    {
        public double spacing;
        public double theoSpacing;
        public double maxIntensity;
        public int charge;
        public bool light;
        public bool heavy;
        public double MZ;

        public Spacing(MZPeak light, MZPeak heavy, int charge)
        {
            spacing = (heavy.MZ - light.MZ) * (double)charge;

            if (heavy.Intensity >= light.Intensity)
            {
                maxIntensity = heavy.Intensity;
                this.heavy = true;
            }
            else
            {
                maxIntensity = light.Intensity;
                this.light = true;
            }
        }

        public Spacing(double exp, double theo, int z, double light, double heavy)
        {
            spacing = exp;
            theoSpacing = theo;
            charge = z;
            MZ = (light + heavy) * 0.5;
        }

        public Spacing()
        {
            spacing = double.NaN;
            maxIntensity = double.NaN;
        }

        public int CompareTo(Spacing other)
        {
            int comp = this.spacing.CompareTo(other.spacing);
            return comp;
        }
    }
}
