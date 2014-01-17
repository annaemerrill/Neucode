using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coon.NeuQuant
{
   class PeakPattern
    {
        public int NumPeaks;
        public int PatternID;
        public double RawScore
        {
            get
            {
                double score = 0;
                int peakCount = 0;
                if (PatternPairs != null)
                {
                    foreach (IsotopePair isotopePair in PatternPairs)
                    {
                        for (int i = 0; i < ParentPeptide.numChannels; i++)
                        {
                            if (isotopePair.ChannelPeaks[i] != null)
                            {
                                peakCount++;
                                score += isotopePair.ChannelPeaks[i].GetSignalToNoise() * (2 - isotopePair.Isotope / 2);
                            }
                        }
                    }
                    score = score / peakCount;
                }
                return score;
            }
        }
        public double NormalizedScore { set; get; }
        public PeptideID ParentPeptide;
        public List<IsotopePair> PatternPairs { set; get; }
        int PairCount
        {
            get
            {
                if (PatternPairs != null) return PatternPairs.Count;
                else return 0;
            }
        }

        public PeakPattern(int numPeaks, int patternID, PeptideID parentPeptide, IsotopePair isotopePair = null)
        {
            NumPeaks = numPeaks;
            PatternID = patternID;
            ParentPeptide = parentPeptide;

            PatternPairs = new List<IsotopePair>();

            if (isotopePair != null)
            {
                PatternPairs.Add(isotopePair);
            }

        }

        public void setNormalizedScore(double normalizationValue)
        {
            NormalizedScore = PairCount + (RawScore / normalizationValue);
        }

    }
}
