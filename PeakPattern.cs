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
                    foreach (Pair pair in PatternPairs)
                    {
                        for (int i = 0; i < ParentPeptide.numChannels; i++)
                        {
                            for (int j = 0; j < ParentPeptide.numIsotopes; j++)
                            {
                                if (pair.peaks[i, j] != null)
                                {
                                    peakCount++;
                                    score += pair.peaks[i, j].GetSignalToNoise() * (2 - j / 2);
                                }
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
        public List<Pair> PatternPairs { set; get; }
        int PairCount
        {
            get
            {
                if (PatternPairs != null) return PatternPairs.Count;
                else return 0;
            }
        }

        public PeakPattern(int numPeaks, int patternID, PeptideID parentPeptide, Pair pair = null)
        {
            NumPeaks = numPeaks;
            PatternID = patternID;
            ParentPeptide = parentPeptide;

            PatternPairs = new List<Pair>();

            if (pair != null)
            {
                PatternPairs.Add(pair);
            }

        }

        public void setNormalizedScore(double normalizationValue)
        {
            NormalizedScore = PairCount + (RawScore / normalizationValue);
        }

    }
}
