using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coon.NeuQuant
{
    
    [Flags]
    public enum NonQuantifiableType
    {
        Unspecified = 0,
        NoLabel = 1,
        NotResolvable = 2,
        NoPeaksFoundWithinTolerance = 4,
        WrongIsotopeDistribution = 8,
        ImproperSpacing = 16,
        InconsistentPeakPatterns = 32,
        WrongElutionProfiles = 64,
        NotEnoughMeasurements = 128,
        Quantified = 256,
    }
}
