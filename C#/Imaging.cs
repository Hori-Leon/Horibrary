using Emgu.CV;
using Emgu.CV.Stitching;
using Emgu.CV.Structure;
using Emgu.CV.Util;
using System;
using System.Drawing;
using System.IO;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace Horibrary
{
    class Imaging
    {
        public static RenderTargetBitmap CaptureScreen(BitmapSource RB, ImageSource ReconImage, int[] RectSize_1, int[] RectSize_2, bool RectEn)
        {
            return (RenderTargetBitmap)FrameDrawing(RB, ReconImage, RectSize_1, RectSize_2, 0, RectEn);
        }

        private static object FrameDrawing(BitmapSource RB, ImageSource ReconImage, int[] RectSize_1, int[] RectSize_2, short FeatureFlag, bool RectEn)
        {
            BitmapEncoder encoder = new BmpBitmapEncoder();
            MemoryStream stream = new MemoryStream();

            int C_X = (int)(RB.Width - RectSize_1[0]) / 2 + 5;
            int C_Y = (int)(RB.Height - RectSize_1[1] - RectSize_1[2]) / 2 - 5;

            DrawingVisual DV = new DrawingVisual();
            using (DrawingContext DC = DV.RenderOpen())
            {
                if (RB != null) DC.DrawImage(RB, new System.Windows.Rect(0, 0, RB.Width, RB.Height));
                if (RectEn)
                {
                    DC.DrawRectangle(null, new System.Windows.Media.Pen(System.Windows.Media.Brushes.Crimson, 2), new System.Windows.Rect(C_X, C_Y, RectSize_1[0], RectSize_1[1]));
                    DC.DrawRectangle(null, new System.Windows.Media.Pen(System.Windows.Media.Brushes.Yellow, 2), new System.Windows.Rect(C_X - 2, C_Y - 2, RectSize_2[0] + 1, RectSize_2[1] + 1));
                }
                if (ReconImage != null) DC.DrawImage((BitmapSource)ReconImage,
                    new System.Windows.Rect(C_X, C_Y, RectSize_1[0], RectSize_1[1]));
            }
            RenderTargetBitmap RTBM = new RenderTargetBitmap((int)RB.Width, (int)RB.Height, 96, 96, PixelFormats.Default);
            RTBM.Render(DV);

            if (FeatureFlag == 1)
            {
                encoder.Frames.Add(BitmapFrame.Create(RTBM));
                encoder.Save(stream);
                return stream;
            }
            else
            {
                stream.Flush();
                return RTBM;
            }
        }

        public static BitmapImage Bitmap2BitampImg(Bitmap src)
        {
            MemoryStream ms = new MemoryStream();
            ((Bitmap)src).Save(ms, System.Drawing.Imaging.ImageFormat.Bmp);
            BitmapImage image = new BitmapImage();
            image.BeginInit();
            ms.Seek(0, SeekOrigin.Begin);
            image.StreamSource = ms;
            image.EndInit();
            return image;
        }

        public static BitmapImage Byte2BitampImg(byte[] Arr)
        {
            using (var Stream = new System.IO.MemoryStream(Arr))
            {
                var image = new BitmapImage();
                image.BeginInit();
                image.CacheOption = BitmapCacheOption.OnLoad;
                image.StreamSource = Stream;
                image.EndInit();
                return image;
            }
        }


        private void Stitching(ImageSource[] IMGs)
        {
            VectorOfMat VectM = new VectorOfMat();

            int Xoffset = 80;
            int Yoffset = 80;

            for(int i = 0; i < IMGs.Length; i++)
            {
                using (var Pic = new MemoryStream())
                {
                    BitmapEncoder Enco = new BmpBitmapEncoder();
                    byte[] StreamBytes = null;

                    Enco.Frames.Add((BitmapFrame)IMGs[i]);
                    Enco.Save(Pic);
                    StreamBytes = Pic.ToArray();

                    using (Mat Back2Mat = new Mat()) // Frame image buffer
                    {
                        CvInvoke.Imdecode(StreamBytes, Emgu.CV.CvEnum.ImreadModes.AnyColor, Back2Mat);
                        VectM.Push(Back2Mat);
                    }

                    Pic.Flush();
                }
            }

            using (Stitcher ST = new Stitcher(Stitcher.Mode.Panorama))
            {
                Mat result = new Mat();
                ST.Stitch(VectM, result);
                byte[] ImgStream = result.ToImage<Bgr, byte>().ToJpegData();
                CroppedBitmap CB = new CroppedBitmap(Byte2BitampImg(ImgStream), new Int32Rect(Xoffset, Yoffset, result.Width - Xoffset * 2, result.Height - Yoffset * 2));

                using (FileStream FS = File.OpenWrite("Panorama_Result.png"))
                {
                    BitmapEncoder encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(CB));
                    encoder.Save(FS);
                    FS.Flush();

                    string Path = System.AppDomain.CurrentDomain.BaseDirectory;
                }
                ST.Dispose();
                result.Dispose();
            }
            GC.Collect();
        }

    }
}
