using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OMNE
{
    public class PSM : IComparable<PSM>
    {

        public string RawFileName;

        public int ScanNumber;

        public PSM(string rawfileName, int scanNumber)
        {
            RawFileName = rawfileName;
            ScanNumber = scanNumber;
        }        

        public int CompareTo(PSM other)
        {
            int comp = this.RawFileName.CompareTo(other.RawFileName);
            if (comp == 0)
            {
                return ScanNumber.CompareTo(other.ScanNumber);
            }
            else
            {
                return comp;
            }        
        }
    }
}
