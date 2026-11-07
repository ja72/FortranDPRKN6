using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA
{
    public struct Vector : 
        ICollection<double>, 
        System.Collections.ICollection, 
        IEquatable<Vector> 
    {
        readonly Memory<double> data;

        #region Factory
        public Vector(int size)
        {
            data = new double[size];
        }
        public Vector(int size, Func<int, double> initializer)
        {
            var data = new double[size];
            for (int i = 0; i < size; i++)
                data[i] = initializer(i);
            this.data = new Memory<double>(data);
        }
        public Vector(ReadOnlySpan<double> data)
        {
            this.data = new Memory<double>(data.ToArray());
        }
        public Vector(Memory<double> data)
        {
            this.data = data;
        }
        public Vector(params double[] data)
        {
            this.data = new Memory<double>(data);
        } 
        #endregion

        #region Properties
        public ref double this[int index]
        {
            get => ref data.Span[index];
        }
        public Vector this[int index, int length]
        {
            get => new Vector(data.Slice(index, length));
        }
        public int Length => data.Length;

        public double SumSquares => Dot(this, this);
        public double Magnitude => Math.Sqrt(SumSquares); 
        #endregion

        #region Algebra
        public static Vector Negate(Vector a)
            => new Vector(
                a.Length, 
                (idx) => -a[idx]);
        public static Vector Scale(double factor, Vector a)
            => new Vector(
                a.Length, 
                (idx) => factor * a[idx]);
        public static Vector Add(Vector a, Vector b)
            => new Vector(
                a.Length, 
                (idx) => a[idx] + b[idx]);
        public static Vector Subtract(Vector a, Vector b)
            => new Vector(
                a.Length, 
                (idx) => a[idx] - b[idx]);

        public static Vector operator +(Vector a, Vector b) => Add(a, b);
        public static Vector operator -(Vector a) => Negate(a);
        public static Vector operator -(Vector a, Vector b) => Subtract(a, b);
        public static Vector operator *(double f, Vector a) => Scale(f, a);
        public static Vector operator *(Vector a, double f) => Scale(f, a);
        public static Vector operator /(Vector a, double d) => Scale(1/d, a);
        #endregion

        #region Linear Algebra

        public static double Dot(Vector a, Vector b)
        {
            int n = Math.Min(a.Length, b.Length);
            double sum = 0.0;
            for (int i = 0; i < n; i++)
                sum += a[i] * b[i];
            return sum;
        }

        public static Vector Normalize(Vector a)
        {
            double mag = a.Magnitude;
            if (mag == 0.0)
                throw new InvalidOperationException("Cannot normalize a zero vector.");
            return Scale(1.0 / mag, a);
        }

        #endregion

        #region ICollection
        public bool IsReadOnly
        {
            get => true;
        }
        public int Count
        {
            get => data.Length;
        }
        public double[] ToArray()
        {
            return data.ToArray();
        }
        public int IndexOf(double item) => throw new NotSupportedException();
        public bool Contains(double item) => IndexOf(item)>=0;
        public void CopyTo(Array array, int index)
            => Array.Copy(ToArray(), 0, array, index, Count);
        public void CopyTo(double[] array, int arrayIndex)
            => CopyTo(array as Array, arrayIndex);
        public IEnumerator<double> GetEnumerator()
        {
            var span = data.ToArray();
            for (int i = 0; i < data.Length; i++)
                yield return span[i];   
        }
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
        bool System.Collections.ICollection.IsSynchronized
        {
            get => false;
        }
        object System.Collections.ICollection.SyncRoot
        {
            get => new NotSupportedException();
        }
        void ICollection<double>.Add(double item) => throw new NotSupportedException();
        void ICollection<double>.Clear() => throw new NotSupportedException();
        bool ICollection<double>.Remove(double item) => throw new NotSupportedException();

        #endregion

        #region Equality
        public bool Equals(Vector other) => data.Span.SequenceEqual(other.data.Span);
        public override bool Equals(object obj) => obj is Vector other && Equals(other);
        public override int GetHashCode() => data.GetHashCode();

        #endregion
    }
}
