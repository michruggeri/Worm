#ifndef RANDOM__HPP
#define RANDOM__HPP

extern double PG;

class Random {
    public:
        Random(long seed=1);
        Random(const Random &model);
        ~Random(){};

        long Seed() const {return _seed;};

        void SetSeed(long newseed);

        int Int_Random(int n, int m) const;
        
        double Double_Random(double n, double m) const;
        double Gauss_Random(double av, double sigma) const;

        const Random &operator=(const Random &random);

    private:
        long _seed;
};

#endif 
