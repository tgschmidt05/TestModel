using System;

namespace ClassLibrary1
{
    public class Class1
    {
        public double addingStuff(double a, double b)
        {
            return a + b;
        }

        public double subtractingStuff(double a, double b)
        {
            return a - b;
        }
    }

    public class Class2
    {
        static void Main(string[] args)
        {
            Class1 test = new Class1();
            double output1 = test.addingStuff(4, 7);
            double output2 = test.subtractingStuff(4, 7);
            
            Console.WriteLine(output1);
            Console.WriteLine(output2);
        }
    }
}