import java.util.Random;

import static java.lang.StrictMath.*;


public class Helper {

    public double exponential_distribution(double lambda){
        Random r = new Random();
        float j1 = (float) r.nextInt(1000+1)/1000;
        if (j1 == 0) {
            j1 += 0.0001;
        }
        float j2 = (float) (-log(j1)/lambda);
        return j2;
    }


    public int poisson_distribution(double lambda){
        float k, L, j1, j2, j3;
        int p;
        Random r = null;
        j1 = (float) r.nextInt(1000+1)/1000;
        k = 0;
        L = (float) exp(-lambda);
        j3 = 0;
        do{
            j2 = (float) (L * pow(lambda, k));
            p = 1;
            for (int i6 = 0; i6 <= k; i6++){
                if (i6 == 0)
                    p = 1;
                else
                    p *= i6;
            }
            j2 /= p;
            j3 += j2;
            k++;
        } while (j1 >= j3);

        return (int) (k-1);
    }

    public int normal_distribution(double mean, double stdev){
        // TO MODEL BASED ON CUMULATIVE DENSITY FUNCTION OF NORMAL DISTRIBUTION BASED ON BOOK OF SHELDON ROSS, Simulation, The polar method, p80.

        float v1, v2, t;
        int x;
        do{
            Random r = new Random();
            v1 = (float) r.nextInt(1000+1)*2;
            v1 /= 1000;
            v1 -= 1;
            v2 = (float) r.nextInt(1000+1)*2;
            v2 /= 1000;
            v2 -= 1;
            t=v1*v1+v2*v2;
        }
        while(t>=1||t==0);
        float multiplier = (float) sqrt(-2*log(t)/t);
        x = (int) (v1 * multiplier * stdev + mean);
        return x;
    }

    public boolean bernouilli_distribution(double prob){
        Random r = new Random();
        float j1 = (float) r.nextInt(1000+1)/1000;
        if (j1 < prob)
            return false;
        else
            return true;
    }

    public int uniform_distribution(double a, double b){
        Random r = new Random();
        float j1 = (float) r.nextInt(1000+1)/1000;
        return (int) ((int) a + (b-a) * j1);
    }


    public int triangular_distribution(int a, int b, int c){
        float mean, stdev,x, L;

        mean = (a+b+c)/3;
        stdev = (float) ((pow(a,2)+pow(b,2)+pow(c,2)-a*b-a*c-b*c)/18);
        stdev = (float) sqrt(stdev);
        Random r = null;
        float j1 = (float) r.nextInt(1000+1)/1000;
        x = a;

        do{
            if (x <= b)
                L = (float) (pow((x-a),2)/((c-a)*(b-a)));
            else
                L = (float) (1-(pow(c-x,2)/((c-a)*(c-b))));
            x++;
        } while (j1 >= L);

        return (int) (x-1);
    }
}