using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon;

namespace HILAC
{
    class NonQuantifiablePeptideID
    {
        public int scanNumber;
        public Peptide sequence;
        public int reasonNQ;

        public NonQuantifiablePeptideID(int scanNumber, Peptide sequence)
        {
            this.scanNumber = scanNumber;
            this.sequence = new Peptide(sequence);
        }

        public NonQuantifiablePeptideID(int scanNumber, Peptide sequence, int reason)
        {
            this.scanNumber = scanNumber;
            this.sequence = new Peptide(sequence);
            reasonNQ = reason;
        }

        public string getReasonNQ()
        {
            if (reasonNQ == 0)
            {
                return "# lysines: 0";
            }
            else if (reasonNQ == 1)
            {
                return "# quantifiable pairs: 0";
            }
            else if (reasonNQ == 2)
            {
                return "not enough quantifiable pairs";
            }
            else if (reasonNQ == 3)
            {
                return "not able to align pairs";
            }
            else if (reasonNQ == 4)
            {
                return "not able to intensity-filter aligned pairs";
            }
            else
            {
                return "not able to quantify intensity-filtered pairs";
            }
        }
    
    }

    

}
