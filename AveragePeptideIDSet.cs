using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HILAC
{
    class AveragePeptideIDSet
    {
        public List<AveragePeptideID> peptides;
        public int totalHILACPairs
        {
            get
            {
                if (peptides != null)
                {
                    int count = 0;
                    for (int i = 0; i < peptides.Count(); i++)
                    {
                        count += peptides[i].hILACPairs.Count;
                    }
                    return count;
                }
                return 0;
            }
        }

        public List<HILACPair> allHILACPairs;
        public List<HILACPair> alignedHILACPairs;
        public List<HILACPair> filteredHILACPairs;
        public List<HILACPair> quantifiedHILACPairs;

        public AveragePeptideIDSet()
        {
            peptides = new List<AveragePeptideID>();
            allHILACPairs = new List<HILACPair>();
            alignedHILACPairs = new List<HILACPair>();
            filteredHILACPairs = new List<HILACPair>();
            quantifiedHILACPairs = new List<HILACPair>();
            
        }

    }
}
