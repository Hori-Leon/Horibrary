package com.gcam.common;

import android.content.Context;
import android.graphics.Bitmap;

import com.github.mikephil.charting.data.Entry;
import com.github.mikephil.charting.data.LineData;
import com.github.mikephil.charting.data.LineDataSet;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;

public class Formula {

    private double[][] Grid_S;
    private double[][] Grid_M;
    private double[][] ER_CS;
    private double[][] ER_BA;
    private double[][] ER_CO;

    int ROW_PIXEL_SIZE = 12;
    int NUCLIDE_NUMBER = 3;
    int CUTTING_MARGIN = 85;
    double MaskRF = 1.24;
    double ImageRF = 3;
    int SizeX = 512;
    int SizeY = 512;

    /**
     * Test Variables
     */
    private double[][] TEST_ADC_M;
    private double[][] TEST_ADC_A;
    public int[] TEST_RESULT;
    private boolean[] TEST_FT_BOOL;
    private int[] TEST_RESULT_T;
    public Bitmap TEST_IMAGE;

    public Formula(Context mContext)
    {
        /**
         * Initialize Label, Mask Pattern, Nuclides' Windows
         */
        Grid_S = DataLoad("LABEL.txt",mContext,"\t");
        Grid_M = DataLoad("19x19.txt",mContext,"\t");
        ER_CS = DataLoad("ER_CS.txt",mContext,"\t");
        ER_BA = DataLoad("ER_BA.txt",mContext,"\t");

        ER_CO = new double[2][144];
        for(int i = 0; i< 2; i++){
            for(int j = 0; j<144; j++){
                ER_CO[i][j] = ER_CS[i][j] * 1100 / 662;
            }
        }
        Grid_M = Bi_Inter(Grid_M, MaskRF);

/*
        // Test Area.
        TEST_ADC_M = LoadADC("M.txt",mContext);
        TEST_ADC_A = LoadADC("A.txt",mContext);
        TEST_FT_BOOL = new boolean[5];
        TEST_FT_BOOL[0] = true;
        TEST_RESULT = ReconGammaImage(TEST_ADC_M,TEST_ADC_A,TEST_FT_BOOL);
        int a = (int) sqrt(TEST_RESULT.length);
        TEST_IMAGE = GetPixelARGB(TEST_RESULT,200,10);
        */
    }

    /**
     * Load Initializing data from asset.
     * @param Filename
     * @return
     */
    private double[][] DataLoad(String Filename, Context mContext, String Regex) {
        BufferedReader reader = null;

        String mLine;
        String[] Row;
        int TOTAL_LENGTH = 0;
        int COL= 0;
        int ROW = 0;

        double[] Buf = null;
        double[][] Out = null;

        try {
            reader = new BufferedReader(
                    new InputStreamReader(mContext.getAssets().open(Filename), "UTF-8"));


            while ((mLine = reader.readLine()) != null) {
                Row = mLine.split(Regex);

                if (Buf == null)
                    if(Regex == "\t") Buf = new double[Row.length * 1000]; // Initialize Buf Approximately. Because i don't know the size of text file.
                else Buf = new double[Row.length * 50000];

                for (int i = 0; i < Row.length; i++)
                    Buf[i + Row.length * COL] = Double.parseDouble(Row[i]);


                COL++;
                TOTAL_LENGTH += Row.length;
            }
            ROW = TOTAL_LENGTH / COL;

        } catch (IOException e) {
            //log the exception
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                    for (int i = 0; i <  ROW; i++) {
                        for (int j = 0; j < COL; j++) {
                            if((ROW > COL) || (ROW == COL)){
                                if(Out ==  null)Out = new double[COL][ROW];
                                Out[j][i] = Buf[j * ROW + i];
                            }
                            else{
                                if(Out ==  null) Out = new double[ROW][COL];
                                Out[i][j] = Buf[i * COL + j];
                            }
                        }
                    }
                } catch (IOException e) {
                    //log the exception
                }
            }
        }

        return Out;
    }

    /**
     * Load ADC Data for debugging. Can be ignored or deprecated.
     * @param Filename
     * @param mContext
     * @return
     */
    private double[][] LoadADC(String Filename, Context mContext){

        BufferedReader reader = null;

        String mLine;
        String[] Row;
        int TOTAL_LENGTH = 0;
        int TOTAL_LINE = 0;

        double[] Buf = null;
        double[][] Out = null;

        try {
            reader = new BufferedReader(
                    new InputStreamReader(mContext.getAssets().open(Filename), "UTF-8"));


            while ((mLine = reader.readLine()) != null) {
                Row = mLine.split(" ");

                if (Buf == null) Buf = new double[Row.length * 50000];

                for (int i = 0; i < Row.length; i++)
                    Buf[i + Row.length * TOTAL_LINE] = Integer.parseInt(Row[i]);


                TOTAL_LINE++;
                TOTAL_LENGTH += Row.length;
            }

        } catch (IOException e) {
            //log the exception
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                    Out = new double[4][TOTAL_LINE];
                    for (int i = 0; i < 4 ; i++) {
                        for (int j = 0; j < TOTAL_LINE ; j++) {
                            Out[i][j] = Buf[(j * 4) + i];
                        }
                    }
                } catch (IOException e) {
                    //log the exception
                }
            }
        }



        return Out;
    }

    /**
     * Second Dimensional Cross Correlation. Mainly used.
     * @param Dim1
     * @param Dim2
     * @return
     */
    public double[][] xcorr2(double[][] Dim1, double[][] Dim2)
    {
        int Dim1_H = Dim1[0].length;
        int Dim1_W = Dim1[1].length;
        int Dim2_H = Dim2[0].length;
        int Dim2_W = Dim2[1].length;

        // Xt -> Zero - Padding
        double[][] Xt = new double[Dim1_H * 2 + Dim2_H - 2][Dim1_W * 2 + Dim2_W - 2];
        double[][] Out = new double[Xt[0].length - Dim1_H + 1][Xt[1].length - Dim1_W + 1];
        double Foo = 0;

        // Fill Dim1 Values into the center of Xt.
        for (int i = Dim1_H - 1; i < Dim1_H + Dim2_H - 1; i++)
        {
            for (int j = Dim1_W - 1; j < Dim2_W + Dim1_W - 1; j++)
            {
                Xt[i][j] = Dim2[i - (Dim1_H - 1)][j - (Dim1_W - 1)];
            }
        }

        // Multiply every values in both Dim1 & Dim2 while moving Dim2
        for (int i = 0; i < Out[0].length; i++)
        {
            for (int j = 0; j < Out[1].length; j++)
            {
                for (int k = 0; k < Dim1_H; k++)
                {
                    for (int l = 0; l < Dim1_W; l++)
                    {
                        Foo += Dim1[k][l] * Xt[i + k][j + l];
                    }
                }
                Out[i][j] = Foo;
                Foo = 0;
            }
        }

        return Out;

    }

    /**
     *  Second Dimensional Convolution. Deprecated but might be used in future. Same result as xcorr2
     * @param Dim1
     * @param Dim2
     * @return
     */
    public double[][] conv2(double[][] Dim1, double[][] Dim2)
    {
        int Dim1_H = Dim1[0].length;
        int Dim1_W = Dim1[1].length;
        int Dim2_H = Dim2[0].length;
        int Dim2_W = Dim2[1].length;

        int Center_H = (int) floor((double)(Dim2_H + 1) / 2);
        int Center_W = (int) floor((double)(Dim2_W + 1) / 2);

        // Four Location Pixels of Center
        int Center_L = Center_W - 1;
        int Center_R = Dim2_W - Center_W;
        int Center_U = Center_H - 1;
        int Center_D = Dim2_H - Center_H;

        double[][] Ref = new double[Dim1_H + Center_U + Center_D][Dim1_W + Center_L + Center_R];

        for (int X = Center_U; X < Dim1_H + Center_U; X++)
        {
            for (int Y = Center_L; Y < Center_L + Dim1_W; Y++)
            {
                Ref[X][Y] = Dim1[X - Center_U][Y - Center_L];
            }
        }

        double[][] Out = new double[Dim1_H][Dim1_W];
        for (int X = 0; X < Dim1_H; X++)
        {
            for (int Y = 0; Y < Dim1_W; Y++)
            {
                for (int i = 0; i < Dim2_H; i++)
                {
                    for (int j = 0; j < Dim2_W; j++)
                    {
                        Out[X][Y] += (Ref[i + X][j + Y] * Dim2[i][j]);
                    }
                }
            }
        }

        return Out;
    }

    /**
     * Bilinear Interpolation from MATLAB, Porting in C# Simply. (imresize)
     * Although it does not provide same results as MATLAB, Similarity rate of Result is about 99% above.
     * @param mMat
     * @param Scale
     * @return
     */
    private double[][] Bi_Inter(double[][] mMat, double Scale)
    {
        int Height = mMat.length;
        int Width = mMat.length;

        Height = (int)Math.round(Height * Scale);
        Width = (int)Math.round(Width * Scale);

        double[][] Buf = new double[Height][Width];

        for (int i = 1; i < Height + 1; i++)
        {
            double Y = (1 / Scale * i) + (0.5 * (1 - 1 / Scale));
            for (int j = 1; j < Width + 1; j++)
            {
                double X = (1 / Scale * j) + (0.5 * (1 - 1 / Scale));
                if (X < 1) X = 1;
                else if (X > mMat.length - 0.001) X = mMat.length - 0.001;
                if (Y < 1) Y = 1;
                else if (Y > mMat.length - 0.001) Y = mMat.length - 0.001;
                int X_1 = (int)Math.floor(X);
                int Y_1 = (int)Math.floor(Y);
                int X_2 = X_1 + 1;
                int Y_2 = Y_1 + 1;

                // Calculate each location of 4 Neighbor Pixels
                double NP_1 = mMat[Y_1 - 1][X_1 - 1];
                double NP_2 = mMat[Y_1 - 1][X_2 - 1];
                double NP_3 = mMat[Y_2 - 1][X_1 - 1];
                double NP_4 = mMat[Y_2 - 1][X_2 - 1];

                // Calculate each weight of 4 Neighbor Pixels
                double PW_1 = (Y_2 - Y) * (X_2 - X);
                double PW_2 = (Y_2 - Y) * (X - X_1);
                double PW_3 = (X_2 - X) * (Y - Y_1);
                double PW_4 = (Y - Y_1) * (X - X_1);

                Buf[i - 1][j - 1] = PW_1 * NP_1 + PW_2 * NP_2 + PW_3 * NP_3 + PW_4 * NP_4;
            }
        }

        return Buf;
    }

    /**
     * Rotate Matrix for 90 degree
     * @param mMat
     * @return
     */
    private double[][] Rot90(double[][] mMat)
    {
        int Mat_H = mMat[0].length;
        int Mat_W = mMat[1].length;
        double[][] Out = new double[Mat_H][Mat_W];

        for (int X = 0; X < Mat_H; X++)
        {
            int X_2 = Mat_H - 1 - X;
            for (int Y = 0; Y < Mat_W; Y++)
            {
                int Y_2 = Mat_W - 1 - Y;
                Out[X][Y] = mMat[X_2][Y_2];
            }
        }
        return Out;
    }

    /**
     * Rotate Matrix for 90 degree CounterClockwise
     * @param mMat
     * @return
     */
    private double[][] Rot90CC(double[][] mMat)
    {
        int Mat_H = mMat[0].length;
        int Mat_W = mMat[1].length;
        double[][] Out = new double[Mat_H][Mat_W];

        for (int X = 0; X < Mat_H; X++)
        {
            for (int Y = 0; Y < Mat_W; Y++)
            {
                int Y_2 = Mat_W - 1 - Y;
                Out[X][Y] = mMat[Y_2][X];
            }
        }
        return Out;
    }


    /**
     * Transpose Matrix
     * @param Input
     * @return
     */
    private double[][] Transpose(double[][] Input)
    {
        double[][] Out = new double[Input[0].length][Input[1].length];

        for (int i = 0; i < Input[0].length; i++)
        {
            for (int j = 0; j < Input[1].length; j++)
            {
                Out[j][i] = Input[i][j];
            }
        }

        return Out;
    }

    /**
     * Flip Maxtrix from Left to Right.
     * @param Input
     * @return
     */
    private double[][] Flip_LR(double[][] Input)
    {
        double[][] Out = new double[Input[0].length][Input[1].length];

        for (int i = 0; i < Input[0].length; i++)
        {
            for (int j = 0; j < Input[1].length; j++)
            {
                Out[i][Input[1].length - j - 1] = Input[i][j];
            }
        }
        return Out;
    }

    /**
     * Histogram Construction. Check commentout before using first.
     * @param Data
     * @param Bin_size
     * @param Bin_Limit
     * @return
     */
    public LineData Histogram(int[] Data, int Bin_size, int Bin_Limit)
    {
        // Bin_Limit must be larger than Bin_Size.
        // Bin_Limit must be divided by Bin_Size totally in Static-Ranged Histogram Mode.
        // Bin_Limit can be set automatically in Dynamic-Ranged Histogram Mode.
        // Please remind this before using.

        // If you want to use Dynamic-Ranged Histogram, set Value to 1. Otherwise, 0;
        int Value = 0;
        int Parity = 0;



        // Static-Ranged Histogram
        if (Bin_Limit == 0 || Bin_Limit % Bin_size != 0)
        {
            if (Value == 0) return null;
                // Dynamic-Ranged Histogram.
            else if (Value == 1)
            {
                Arrays.sort(Data);
                Bin_Limit = ((Data[Data.length - 1] / Bin_size) + 1) * Bin_size;
            }
        }

        double[] Out = new double[Bin_Limit / Bin_size * 2];

        for (int i = 0; i < Bin_Limit / Bin_size; i++)
        {
            Out[i] = i * Bin_size;
        }

        for (int i = 0; i < Data.length; i++)
        {
            if ((int)Data[i] >= Bin_Limit) Out[Out.length - 1]++;
            else Out[(Out.length / 2) + Data[i] / Bin_size]++;
        }

        // Fill Zero-Padding
        Out[Out.length - 1] = 0;
        Out[0] = 0;

        // LineChart Integration
        return MakeLineChart(Out,"Histo");
    }

    public LineData MakeLineChart(double[] InputArr,String Datatype){

        // Output
        ArrayList<Entry> Values = new ArrayList<>();
        LineDataSet OutputLineSet;
        LineData OutputLine;

        if(Datatype.equals("Histo")){
            for(int i = 0; i < InputArr.length / 2; i++){
                Values.add(new Entry( (float)InputArr[i], (float) InputArr[InputArr.length / 2 + i] ));
            }
        } else if(Datatype.equals("Normal")){
            for(int i = 0; i < InputArr.length / 2; i++){
                Values.add(new Entry( i, (float)InputArr[i]) );
            }
        }

        OutputLineSet = new LineDataSet(Values,null);
        OutputLineSet.setLineWidth(2f);
        OutputLineSet.setDrawCircles(false);
        OutputLineSet.setHighlightEnabled(false);
        OutputLine = new LineData(OutputLineSet);
        OutputLine.setDrawValues(false);

        return OutputLine;
    }

    /**
     * Stack-up ADC values.
     * @param A
     * @param B
     * @return
     */
    public double[][] ADCDataCombine(double[][] A, double[][] B)
    {
        double[][] Out = new double[4][A[0].length + B[0].length];

        // ADC A Combine
        System.arraycopy(A[0], 0, Out[0], 0, A[0].length);
        System.arraycopy(B[0], 0, Out[0], A[0].length, B[0].length);

        System.arraycopy(A[1], 0, Out[1], 0, A[1].length);
        System.arraycopy(B[1], 0, Out[1], A[1].length, B[1].length);

        System.arraycopy(A[2], 0, Out[2], 0, A[2].length);
        System.arraycopy(B[2], 0, Out[2], A[2].length, B[2].length);

        System.arraycopy(A[3], 0, Out[3], 0, A[3].length);
        System.arraycopy(B[3], 0, Out[3], A[3].length, B[3].length);
        return Out;
    }

    /**
     * Stack-up Energy values.
     * @param A
     * @param B
     * @return
     */
    public int[] EnergyDataCombine(int[] A, int[] B)
    {
        int[] Out = new int[A.length + B.length];

        // Energy Value Combine
        System.arraycopy(A, 0, Out, 0, A.length);
        System.arraycopy(B, 0, Out, A.length, B.length);

        return Out;
    }

    /**
     * Jet Colorscale Constructor. 8-bit ARGB Indexed.
     * @param AlphaIntense
     * @param Fill_Factor
     * @return
     */
    public int[][] JetPalette(int AlphaIntense, int Fill_Factor)
    {
        // 8-bit(256) ARGB indexed Color
        int[][] Out = new int[4][256];

        for (int i = 0; i < Fill_Factor; i++)
        {
            Out[1][255 - i] = 255;

            Out[1][255 - Fill_Factor - i] = 255;
            Out[2][255 - Fill_Factor - i] = 255;

            Out[2][255 - Fill_Factor * 2 - i] = 255;
            Out[3][255 - Fill_Factor * 2 - i] = 255;

            Out[3][255 - Fill_Factor * 3 - i] = 255;
        }

        for (int i = 255; i > 255 - (Fill_Factor * 4); i--) Out[0][i] = AlphaIntense;
        // Debug
        //for (int i = 0; i < 256; i++) Out[0][i] = AlphaIntense;

        return Out;
    }

    /**
     * Debug Palette. Deprecated
     * @return
     */
    public int[][] TestPalette()
    {
        // 8-bit(256) ARGB indexed Color
        int[][] Out = new int[4][256];

        for (int i = 0; i < 256; i++)
        {
            Out[0][i] = 60;
            Out[1][i] = 127;
            Out[2][i] = 17;
        }

        //for (int i = 255; i > 255 - (Fill_Factor * 4); i--) Out[0][i] = AlphaIntense;
        // Debug
        for (int i = 0; i < 256; i++) Out[3][i] = 100;

        return Out;
    }

    /**
     * Conversion from int[] to Bitmap scale with Jet Palette
     * @param PixelData
     * @param AlphaIntense
     * @param Fill_Factor
     * @return
     */
    public Bitmap GetPixelARGB(int[] PixelData, int AlphaIntense, int Fill_Factor) {

        int[] Convert = new int[PixelData.length];
        Bitmap ResultPic = null;
        //int[][] Pal = TestPalette();
        int[][] Pal = JetPalette(AlphaIntense,Fill_Factor);

        for(int i = 0; i < PixelData.length; i++){
            Convert[i] = Pal[0][(int)PixelData[i]] << 24 | Pal[1][(int)PixelData[i]] << 16 |
                    Pal[2][(int)PixelData[i]] << 8 | Pal[3][(int)PixelData[i]] ;
        }

        int Len = (int) sqrt(PixelData.length);

        ResultPic = Bitmap.createBitmap(Convert,Len,Len,Bitmap.Config.ARGB_8888);

        return ResultPic;
    }
}
