using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HILAC
{
    //A light or heavy version of an isotope of a HILAC pair
    class MetabolicForm
    {
        public HILACPairIsotope pair;
        public int scanNumber;
        public double intensity;

        public MetabolicForm(HILACPairIsotope pair, double intensity)
        {
            this.pair = pair;
            this.intensity = intensity;
            scanNumber = pair.pair.scanNumber;
        }
    }
}
