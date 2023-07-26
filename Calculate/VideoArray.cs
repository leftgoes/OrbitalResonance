using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace Calculate
{
    public class DoubleRange
    {
        public readonly double min;
        public readonly double max;
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

        private int LinMap(double x, double f1, double f2, int t1, int t2)
        {
            return (int)((t1 - t2) * (x - f1) / (f2 - f1) + t1);
        }

        public void AddPoint(int frame, double x, double y)
        {
            int i = LinMap(x, xRange.min, xRange.max, 0, width);
            int j = LinMap(y, yRange.min, yRange.max, height, 0);

            if (0 <= i && i < width && 0 <= j && j < height)
            {
                array[frame, j, i]++;
            }
        }

        public byte[] ToByteArray(int frame)
        {
            double arrayMax = array.Cast<double>().Max();
            byte[] bytes = new byte[stride * height];
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    bytes[y * stride + x] = (byte)(255 * array[frame, y, x] / arrayMax);
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
            bmp.Save(filename);
        }
    }
}
