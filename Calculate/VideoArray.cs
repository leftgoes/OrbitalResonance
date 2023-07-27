using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace Calculate
{
    public class DoubleRange
    {
        public readonly double min;
        public readonly double max;

        public DoubleRange()
        {
            this.min = 0;
            this.max = 1;
        }

        public DoubleRange(double min, double max)
        {
            this.min = min;
            this.max = max;
        }
    }

    public class VideoArray
    {
        public readonly int frames;
        public readonly int width;
        public readonly DoubleRange xRange;
        public readonly DoubleRange yRange;
        public readonly int height;
        public double[,,] array;

        public VideoArray(int frames, DoubleRange xRange, DoubleRange yRange, int width, int height) { 
            this.frames = frames;
            this.xRange = xRange;
            this.yRange = yRange;
            this.width = width;
            this.height = height;
            array = new double[frames, height, width];
        }

        public int stride { get { return (width % 4 == 0) ? width : width + 4 - width % 4; } }

        private double LinMap(double x, double f1, double f2, double t1, double t2)
        {
            return (int)((t1 - t2) * (x - f1) / (f2 - f1) + t1);
        }

        private void AddPointInt(int frame, int i, int j, double value)
        {
            if (i < 0 || i >= width || j < 0 || j >= height)
                return;

            Console.WriteLine("Once");
            array[frame, j, i] += value;
        }

        public void AddPoint(int frame, double x, double y, double value = 1)
        {
            double i = LinMap(x, xRange.min, xRange.max, 0, width);
            double j = LinMap(y, yRange.min, yRange.max, height, 0);

            int iInt = (int)i;
            int jInt = (int)j;
            double iFrac = i - iInt;
            double jFrac = j - jInt;
            
            AddPointInt(frame, iInt, jInt, (1 - iFrac) * (1 - jFrac) * value);
            AddPointInt(frame, iInt, jInt + 1, (1 - iFrac) * jFrac * value);
            AddPointInt(frame, iInt + 1, jInt, iFrac * (1 - jFrac) * value);
            AddPointInt(frame, iInt + 1, jInt + 1, iFrac * jFrac * value);
        }

        public byte[] ToByteArray(int frame)
        {
            double arrayMax = array.Cast<double>().Max();
            byte[] bytes = new byte[stride * height];
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    bytes[y * stride + x] = (byte)Math.Round(255 * array[frame, y, x] / arrayMax);
                }
            }

            return bytes;
        }

        public Bitmap ToBitmap(int frame)
        {
            PixelFormat formatOutput = PixelFormat.Format8bppIndexed;
            Rectangle rect = new(0, 0, width, height);

            Bitmap bmp = new(stride, height, formatOutput);
            BitmapData bmpData = bmp.LockBits(rect, ImageLockMode.ReadOnly, formatOutput);

            byte[] bytes = ToByteArray(frame);
            Marshal.Copy(bytes, 0, bmpData.Scan0, bytes.Length);
            bmp.UnlockBits(bmpData);

            return bmp;
        }

        public void SaveFrame(string filename, int frame)
        {
            Bitmap bmp = ToBitmap(frame);
            bmp.Save(filename, ImageFormat.Png);
        }

        public void SaveFrames(string directory)
        {
            if (!Directory.Exists(directory))
                Directory.CreateDirectory(directory);

            for (int frame = 0; frame < frames; frame++)
            {
                SaveFrame(Path.Join(directory, $"frm{frame:05d}.png"), frame);
            }
        }
    }
}
