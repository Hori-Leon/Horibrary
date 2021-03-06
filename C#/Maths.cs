using System;
using System.Collections.Generic;
using System.Linq;


namespace Horibrary
{
    class Maths
    {
        // Bilinear Interpolation Method from MATLAB, Porting in C# Simply. (imresize)
        // Although it does not provide same results as MATLAB, Similarity rate of Result is about 99% above.
        private static double[,] BiInterPol(double[,] Matrix, double Scale)
        {
            int Height = Matrix.GetLength(0);
            int Width = Matrix.GetLength(1);

            Height = (int)Math.Round(Height * Scale);
            Width = (int)Math.Round(Width * Scale);

            double[,] Buf = new double[Height, Width];
            double[,] Out = new double[Height, Width];

            for (int i = 1; i < Height + 1; i++)
            {
                double Y = (1 / Scale * i) + (0.5 * (1 - 1 / Scale));
                for (int j = 1; j < Width + 1; j++)
                {
                    double X = (1 / Scale * j) + (0.5 * (1 - 1 / Scale));
                    if (X < 1) X = 1;
                    else if (X > Math.Sqrt(Matrix.Length) - 0.001) X = Math.Sqrt(Matrix.Length) - 0.001;
                    if (Y < 1) Y = 1;
                    else if (Y > Math.Sqrt(Matrix.Length) - 0.001) Y = Math.Sqrt(Matrix.Length) - 0.001;
                    int X_1 = (int)Math.Floor(X);
                    int Y_1 = (int)Math.Floor(Y);
                    int X_2 = X_1 + 1;
                    int Y_2 = Y_1 + 1;

                    // Calculate each location of 4 Neighbor Pixels
                    double NP_1 = Matrix[Y_1 - 1, X_1 - 1];
                    double NP_2 = Matrix[Y_1 - 1, X_2 - 1];
                    double NP_3 = Matrix[Y_2 - 1, X_1 - 1];
                    double NP_4 = Matrix[Y_2 - 1, X_2 - 1];

                    // Calculate each weight of 4 Neighbor Pixels
                    double PW_1 = (Y_2 - Y) * (X_2 - X);
                    double PW_2 = (Y_2 - Y) * (X - X_1);
                    double PW_3 = (X_2 - X) * (Y - Y_1);
                    double PW_4 = (Y - Y_1) * (X - X_1);

                    Buf[i - 1, j - 1] = PW_1 * NP_1 + PW_2 * NP_2 + PW_3 * NP_3 + PW_4 * NP_4;
                }
            }
            return Out;
        }


        // 2 Dimensional Cross-Correlation Method from MATLAB
        public double[,] Correlation_2D(double[,] Dim1, double[,] Dim2)
        {

            int Dim1_H = Dim1.GetLength(0);
            int Dim1_W = Dim1.GetLength(1);
            int Dim2_H = Dim2.GetLength(0);
            int Dim2_W = Dim2.GetLength(1);

            double[,] Xt = new double[Dim1_H * 2 + Dim2_H - 2, Dim1_W * 2 + Dim2_W - 2];
            for (int i = Dim1_H - 1; i < Dim1_H + Dim2_H - 1; i++)
            {
                for (int j = Dim1_W - 1; j < Dim2_W + Dim1_W - 1; j++)
                {
                    Xt[i, j] = Dim2[i - (Dim1_H - 1), j - (Dim1_W - 1)];
                }
            }

            double[,] Out = new double[Xt.GetLength(0) - Dim1_H + 1, Xt.GetLength(1) - Dim1_W + 1];
            double Foo = new double();
            for (int i = 0; i < Out.GetLength(0); i++)
            {
                for (int j = 0; j < Out.GetLength(1); j++)
                {
                    for (int k = 0; k < Dim1_H; k++)
                    {
                        for (int l = 0; l < Dim1_W; l++)
                        {
                            Foo += Dim1[k, l] * Xt[i + k, j + l];
                        }
                    }
                    Out[i, j] = Foo;
                    Foo = 0;
                }
            }

            return Out;

        }


        // 2 Dimensional Convolution Method from MATLAB, Porting in C# Simply. (conv2 / same)
        // Works Perfectly.
        private double[,] Conv_2D(double[,] Dim1, double[,] Dim2)
        {
            int Dim1_H = Dim1.GetLength(0);
            int Dim1_W = Dim1.GetLength(1);
            int Dim2_H = Dim2.GetLength(0);
            int Dim2_W = Dim2.GetLength(1);

            int Center_H = (int)Math.Floor((double)(Dim2_H + 1) / 2);
            int Center_W = (int)Math.Floor((double)(Dim2_W + 1) / 2);

            // Four Location Pixels of Center
            int Center_L = Center_W - 1;
            int Center_R = Dim2_W - Center_W;
            int Center_U = Center_H - 1;
            int Center_D = Dim2_H - Center_H;

            double[,] Ref = new double[Dim1_H + Center_U + Center_D, Dim1_W + Center_L + Center_R];
            for (int X = Center_U; X < Dim1_H + Center_U; X++)
            {
                for (int Y = Center_L; Y < Center_L + Dim1_W; Y++)
                {
                    Ref[X, Y] = Dim1[X - Center_U, Y - Center_L];
                }
            }

            double[,] Dim2_Rot90 = new double[Dim2_H, Dim2_W];
            Dim2_Rot90 = Rot90(Dim2);

            double[,] Out = new double[Dim1_H, Dim1_W];
            for (int X = 0; X < Dim1_H; X++)
            {
                for (int Y = 0; Y < Dim1_W; Y++)
                {
                    for (int i = 0; i < Dim2_H; i++)
                    {
                        for (int j = 0; j < Dim2_W; j++)
                        {
                            Out[X, Y] += (Ref[i + X, j + Y] * Dim2_Rot90[i, j]);
                        }
                    }
                }
            }

            return Out;
        }

        // Rotate Matrix for 90 degree
        private static T[,] Rot90<T>(T[,] Mat)
        {
            int Mat_H = Mat.GetLength(0);
            int Mat_W = Mat.GetLength(1);
            T[,] Out = new T[Mat_H, Mat_W];

            for (int X = 0; X < Mat_H; X++)
            {
                int X_2 = Mat_H - 1 - X;
                for (int Y = 0; Y < Mat_W; Y++)
                {
                    int Y_2 = Mat_W - 1 - Y;
                    Out[X, Y] = Mat[X_2, Y_2];
                }
            }

            return Out;
        }

        // Transpose Matrix
        private static T[,] Transpose<T>(T[,] Input)
        {
            T[,] Out = new T[Input.GetLength(0), Input.GetLength(1)];

            for (int i = 0; i < Input.GetLength(0); i++)
            {
                for (int j = 0; j < Input.GetLength(1); j++)
                {
                    Out[j, i] = Input[i, j];
                }
            }

            return Out;
        }

        // Flip Maxtrix from Left to Right.
        private static T[,] Flip_LR<T>(T[,] Input)
        {
            T[,] Out = new T[Input.GetLength(0), Input.GetLength(1)];

            for (int i = 0; i < Input.GetLength(0); i++)
            {
                for (int j = 0; j < Input.GetLength(1); j++)
                {
                    Out[i, Input.GetLength(1) - j - 1] = Input[i, j];
                }
            }

            return Out;
        }

        private static double[] Histo(int[] Data, int Bin_size, int Bin_Limit)
        {
            // Bin_Limit must be larger than Bin_Size.
            // Bin_Limit must be divided by Bin_Size totally in Static-Ranged Histogram Mode.
            // Bin_Limit can be set automatically in Dynamic-Ranged Histogram Mode.
            // Please check the comments.

            // If you want to use Dynamic-Ranged Histogram, set Value to 1. Otherwise, 0;
            short Value = 0;

            // Static-Ranged Histogram
            if (Bin_Limit == 0 || Bin_Limit % Bin_size != 0)
            {
                if (Value == 0) return null;
                // Dynamic-Ranged Histogram.
                else if (Value == 1)
                {
                    int Bin_AutoLimit = ((Data.Max() / Bin_size) + 1) * Bin_size;
                    Bin_Limit = Bin_AutoLimit;
                }
            }

            double[] Out = new double[Bin_Limit / Bin_size];

            for (int i = 0; i < Data.Length; i++)
            {
                if ((int)Data[i] >= Bin_Limit) Out[Out.Length - 1]++;
                else Out[Data[i] / Bin_size]++;
            }

            // Fill Zero-Padding
            Out[Out.Length - 1] = 0;
            Out[0] = 0;

            return Out;
        }

        // Jet Scale Palette from MATLAB
        public static int[,] GenJetPalette(int Alpha_Impulse, int Fill_Factor)
        {
            // 8-bit(256) ARGB indexed Color
            int[,] Out = new int[4, 256];

            for (int i = 0; i < Fill_Factor; i++)
            {
                Out[1, 255 - i] = 255;

                Out[1, 255 - Fill_Factor * 1 - i] = 255;
                Out[2, 255 - Fill_Factor * 1 - i] = 255;

                Out[2, 255 - Fill_Factor * 2 - i] = 255;
                Out[3, 255 - Fill_Factor * 2 - i] = 255;

                Out[3, 255 - Fill_Factor * 3 - i] = 255;
            }

            for (int i = 255; i > 255 - (Fill_Factor * 4); i--) Out[0, i] = Alpha_Impulse;
            // Debug
            //for (int i = 0; i < 256; i++) Out[0, i] = Alpha_Impulse;

            return Out;
        }


        public static double[] GaussianFit(double[] Histo_Center, double[] Histo_Count, int[] Fit_Range)
        {
            double[] Out = new double[(Histo_Center.Length + 1)];
            double[] Histo_Center_Ranged = new double[Fit_Range[1] - Fit_Range[0]];
            double[] Histo_Count_Ranged = new double[Fit_Range[1] - Fit_Range[0]];
            Array.Copy(Histo_Center, Fit_Range[0], Histo_Center_Ranged, 0, Histo_Center_Ranged.Length);
            Array.Copy(Histo_Count, Fit_Range[0], Histo_Count_Ranged, 0, Histo_Count_Ranged.Length);

            double Std_Drv = StdDeviation(Histo_Center_Ranged); // Sigma
            int[] Index_Range = new int[Histo_Count_Ranged.Length];
            int Index = 0;
            for (int i = 0; i < Histo_Count_Ranged.Length; i++)
            {
                if ((Histo_Count_Ranged.Max() * 0.85) < Histo_Count_Ranged[i])
                {
                    Index_Range[Index] = i;
                    Index++;
                }
            }

            double Amp = 0;
            for (int i = 0; i < Index; i++)
            {
                Amp += Histo_Count_Ranged[Index_Range[i]];
            }
            Amp = (Amp / Index); // Amplitude

            double Avg = 0;
            for (int i = 0; i < Index; i++)
            {
                Avg += Histo_Center_Ranged[Index_Range[i]];
            }
            Avg = (Avg / Index); // Average

            double[] Gauss = new double[Histo_Center.Length];
            for (int i = 0; i < Histo_Center.Length; i++)
            {
                Gauss[i] = Amp * Math.Exp(-(Math.Pow((i - Avg) / Std_Drv, 2))); // Gaussian Function
            }

            // Confirm whether Curve Fitting is adequate.
            double Gauss_Sum = 0;
            for (int i = Fit_Range[0] * 5; i < Fit_Range[1] * 5; i += 5) Gauss_Sum += Gauss[i];
            if (Histo_Count_Ranged.Sum() * 0.9 > Gauss_Sum) return Out;
            else
            {
                Array.Copy(Gauss, 0, Out, 0, Gauss.Length);
                Out[Gauss.Length] = Std_Drv * Math.Sqrt(4 * Math.Log(2)) / Avg * 100; // Full Width at Half Maximum :: FWHM

                return Out;
            }
        }



        public static double StdDeviation(IEnumerable<double> Values)
        {
            double Out = 0;

            if (Values.Count() > 0)
            {
                double Avg = Values.Average();
                double Sum = Values.Sum(D => Math.Pow(D - Avg, 2));
                Out = Math.Sqrt((Sum) / (Values.Count() - 1));
            }

            return Out;
        }

        public List<object> EnergyResoltion(double[] Histo_Center, double[] Histo_Count)
        {
            List<object> Out = new List<object>();

            double Histo_Size = Histo_Center[1] - Histo_Center[0];
            double SampleCount = Histo_Center.Length;

            // CS-137
            int[] Range = new int[2] { (int)Math.Round(SampleCount / Histo_Size * 0.9), (int)Math.Round(SampleCount / Histo_Size * 1.1) };

            double[] Buf = new double[Histo_Center.Length];
            double Compensator = 0;
            int[] Fit_Range = new int[2];
            while (Buf.Sum() == 0)
            {
                // Fit Range
                int Max_Loca = 0;
                for (int j = Range[0]; j <= Range[1]; j++)
                {
                    if (Histo_Count[j] > Histo_Count[Max_Loca]) Max_Loca = j;
                }

                for (int j = 0; j < Histo_Center.Length; j++)
                {
                    if ((Histo_Center[j] >= Math.Round(Histo_Center[Max_Loca] * (0.99 - Compensator))) && Fit_Range[0] == 0) Fit_Range[0] = j;
                    if ((Histo_Center[j] >= Math.Round(Histo_Center[Max_Loca] * (1.01 + Compensator))) && Fit_Range[1] == 0)
                    {
                        Fit_Range[1] = j;
                        break;
                    }
                }
                Buf = GaussianFit(Histo_Center, Histo_Count, Fit_Range);
                if (Buf.Sum() == 0)
                {
                    Compensator += 0.01;
                    Array.Clear(Fit_Range, 0, Fit_Range.Length);
                }
            }

            double[] FitCurve = new double[Buf.Length - 1];
            Array.Copy(Buf, Buf.Length - 1, FitCurve, 0, FitCurve.Length);
            double FWHM = Buf[Buf.Length];

            Out.Add(FitCurve);
            Out.Add(FWHM);

            return Out;
        }

    }
}
