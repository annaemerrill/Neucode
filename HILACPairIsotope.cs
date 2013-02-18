using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HILAC
{
    /* Light-heavy pair one isotope of a HILAC pair
     * Keeps track of light and heavy metabolic forms
     */   
    class HILACPairIsotope
    {
        public HILACPair pair;
        public MetabolicForm light;
        public MetabolicForm heavy;
        public int scanNumber;
        public int isotope;
        public double ratio
        {
            get
            {
                if (light != null && heavy != null)
                {
                    return heavy.intensity / light.intensity;
                }
                return double.NaN;
            }
        }

        public bool noiseBandCapped;

        public HILACPairIsotope(HILACPair pair, int isotope, double lightIntensity, double heavyIntensity)
        {
            noiseBandCapped = false;
            this.pair = pair;
            this.isotope = isotope;
            scanNumber = pair.scanNumber;
            light = new MetabolicForm(this, lightIntensity);
            heavy = new MetabolicForm(this, heavyIntensity);

            //Check to see if max light should be replaced
            if (pair.parentID.maxLight == null)
            {
                pair.parentID.maxLight = light;
            }
            else
            {
                if (light.intensity > pair.parentID.maxLight.intensity)
                {
                    pair.parentID.maxLight = light;
                }
            }

            //Check to see if max heavy should be replaced
            if (pair.parentID.maxHeavy == null)
            {
                pair.parentID.maxHeavy = heavy;
            }
            else
            {
                if (heavy.intensity > pair.parentID.maxHeavy.intensity)
                {
                    pair.parentID.maxHeavy = heavy;
                }
            }
        }
    }
}
