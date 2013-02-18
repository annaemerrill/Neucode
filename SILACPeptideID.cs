using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon;

namespace HILAC
{
    class SILACPeptideID
    {
        public int idScanNumber;
        public bool idLight;
        public bool idHeavy;
        public Peptide peptideSequence;
        public int charge;
        public double theoreticalMass;
        public double experimentalMass;
        public double theoreticalMZ;
        public double experimentalMZ;
        public double lightTotalIntensity;
        public double heavyTotalIntensity;
        public List<SILACPair> sILACPairs;
        public double retentionTime;
        public double minRT;
        public double maxRT;
        public double averageMassDifference;
        public int lysineCount;


        public SILACPeptideID(int scanNumber, int charge, double theoMass, double expMass, Peptide peptide)
        {
            idScanNumber = scanNumber;
            this.charge = charge;
            theoreticalMass = theoMass;
            experimentalMass = expMass;
            peptideSequence = peptide;
            lightTotalIntensity = 0;
            heavyTotalIntensity = 0;
            sILACPairs = new List<SILACPair>();
            theoreticalMZ = (theoreticalMass + (charge * 1.00727638)) / charge;
            experimentalMZ = (experimentalMass + (charge * 1.00727638)) / charge;
            lysineCount = 0;

            if (peptideSequence.Parent.Mass.MonoisotopicMass == theoreticalMass) //Light version
            {
                idLight = true;
            }
            else //Heavy version
            {
                idHeavy = true;
            }


            string peptid = peptideSequence.Sequence;
            char[] residues = peptid.ToCharArray();

            for (int i = 0; i < residues.Length; i++)
            {
                if (residues[i] == 'K')
                {
                    lysineCount++;
                }
            }

        }

    }
}
