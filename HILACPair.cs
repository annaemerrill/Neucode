using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Coon.Spectra;
using CoonThermo.IO;

namespace HILAC
{
    /* Light-heavy pair from a full scan
     * Keeps track of light and heavy at monoisotope & first 3 isotopes
     */
    class HILACPair
    {
        public AveragePeptideID parentID;
        public HILACPairIsotope[] isotopes;
        public List<HILACPairIsotope> quantIsotopes;
        public int scanNumber;
        public double heavyToLightRatio
        {
            get
            {
                if (totalLightIntensity > 0 && totalLightIntensity > 0)
                {
                    return totalHeavyIntensity / totalLightIntensity;
                }
                return double.NaN;
            }
            set { }
        }

        public double totalLightIntensity
        {
            get
            {
                if (isotopes.Count() > 0)
                {
                    double sum = 0;
                    foreach (HILACPairIsotope isotope in isotopes)
                    {
                        if (isotope != null)
                        {
                            sum += isotope.light.intensity;
                        }
                    }
                    return sum;
                }
                return double.NaN;
            }
            set { }
        }

        public double totalHeavyIntensity
        {
            get
            {
                if (isotopes.Count() > 0)
                {
                    double sum = 0;
                    foreach (HILACPairIsotope isotope in isotopes)
                    {
                        if (isotope != null)
                        {
                            sum += isotope.heavy.intensity;
                        }
                    }
                    return sum;
                }
                return double.NaN;
            }
            set { }
        }

        public bool noiseBandCapped;
        

        public HILACPair(AveragePeptideID parent, int scanNumber)
        {
            parentID = parent;
            this.scanNumber = scanNumber;
            isotopes = new HILACPairIsotope[4];
            totalLightIntensity = 0;
            totalHeavyIntensity = 0;
            quantIsotopes = new List<HILACPairIsotope>();
            noiseBandCapped = false;
        }
        
        public void quantifyHILACPairIsotopesStdDev()
        {
            double ratioAverage;
            double stdDev;
            double[] sqDevFromMean = new double[4];
            double min;
            double max;
            double cumSqDevFromMean = 0;
            int n = 0;
            double ratioSum = 0;
            
            //Calculate average ratio across all isotopes in HILACPair
            foreach (HILACPairIsotope isotope in isotopes)
            {
                if (isotope != null)
                {
                    ratioSum += isotope.ratio;
                    n++;
                }
            }
            ratioAverage = ratioSum / (double)n;

            //Calculate standard deviation across all isotopes in HILACPair
            for (int i = 0; i < isotopes.Count(); i++)
            {
                if (isotopes[i] != null)
                {
                    sqDevFromMean[i] = Math.Pow(Math.Abs(isotopes[i].ratio - ratioAverage), 2);
                    cumSqDevFromMean += sqDevFromMean[i];
                }
            }
            stdDev = Math.Sqrt(cumSqDevFromMean / (double)n);

            //Exclude isotopes that are 1.5StdDevs above or below the mean from quantitation
            min = ratioAverage + (1.5 * stdDev);
            max = ratioAverage + (1.5 * stdDev);
            for (int i = 0; i < isotopes.Count(); i++)
            {
                if (isotopes[i] != null)
                {
                    if (isotopes[i].ratio <= max || isotopes[i].ratio >= min)
                    {
                        quantIsotopes.Add(isotopes[i]);
                        totalLightIntensity += isotopes[i].light.intensity;
                        totalHeavyIntensity += isotopes[i].heavy.intensity;
                    }
                }
            }

            /*
            //Set heavy to light ratio of HILACPair
            if (quantIsotopes.Count == 1)
            {
                heavyToLightRatio = quantIsotopes.ElementAt(0).ratio;
            }
            else if (quantIsotopes.Count == 2)
            {
                heavyToLightRatio = (quantIsotopes.ElementAt(0).ratio + quantIsotopes.ElementAt(1).ratio) / 2.0;
            }
            else if (quantIsotopes.Count == 3)
            {
                heavyToLightRatio = quantIsotopes.ElementAt(1).ratio;
            }
            else //all 4 isotopes quantifiable
            {
                heavyToLightRatio = (quantIsotopes.ElementAt(1).ratio + quantIsotopes.ElementAt(2).ratio) / 2.0;
            }
             */

        }

        //Sets HILAC pair's ratio to the slope between all quantified isotopes
        public void quantifyHILACPairIsotopesLinearRegression()
        {
            for (int i = 0; i < isotopes.Count(); i++)
            {
                if (isotopes[i] != null)
                {
                    quantIsotopes.Add(isotopes[i]);
                    totalLightIntensity += isotopes[i].light.intensity;
                    totalHeavyIntensity += isotopes[i].heavy.intensity;
                }
            }

            if (quantIsotopes.Count == 1)
            {
                heavyToLightRatio = quantIsotopes[0].heavy.intensity / quantIsotopes[0].light.intensity;
            }
            else
            {
                XYPoint[] points = new XYPoint[quantIsotopes.Count];

                for (int i = 0; i < points.Count(); i++)
                {
                    points[i] = new XYPoint(quantIsotopes[i].light.intensity, quantIsotopes[i].heavy.intensity);
                }

                double slope = 0.0;
                double yIntercept = 0.0;
                LeastSquaresFitLinear(points, points.Count(), ref slope, ref yIntercept);

                heavyToLightRatio = slope;
            }
        }

        public static void LeastSquaresFitLinear(XYPoint[] points, int numPoints, ref double M, ref double B)
        {
            //Gives best fit of data to line Y = MC + B  
            double x1, y1, xy, x2, J;
            int i;

            x1 = 0.0;
            y1 = 0.0;
            xy = 0.0;
            x2 = 0.0;

            for (i = 0; i < numPoints; i++)
            {
                x1 = x1 + points[i].X;
                y1 = y1 + points[i].Y;
                xy = xy + points[i].X * points[i].Y;
                x2 = x2 + points[i].X * points[i].X;
            }

            J = ((double)numPoints * x2) - (x1 * x1);
            if (J != 0.0)
            {
                M = (((double)numPoints * xy) - (x1 * y1)) / J;
                M = System.Math.Floor(1.0E3 * M + 0.5) / 1.0E3;
                B = ((y1 * x2) - (x1 * xy)) / J;
                B = System.Math.Floor(1.0E3 * B + 0.5) / 1.0E3;
            }
            else
            {
                M = 0;
                B = 0;
            }
        }  

    }
}
