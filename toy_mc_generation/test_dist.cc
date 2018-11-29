void test_dist() {
auto pT_distr = new TF1("pT_distr", [/*&mass=chi_mass.central*/](double* x, double* p) {

constexpr double beta = 3.39924; // same as in MC generation from Alberto, (was 3.45)
      constexpr double gamma = 0.635858; // same as in MC generation from Alberto, (was 0.73)
      constexpr double A = 1. / (beta - 2.) / gamma;
      const double pT_over_chimass = x[0]/p[0]; // pT = x[0], chimass = p[0]

      return x[0] * pow(1. + A * pT_over_chimass*pT_over_chimass, -beta);
    }, 0,10, 1);

   pT_distr->SetParameter(0, 1);
    pT_distr->Draw();
    std::cout << "max(1) = " << pT_distr->GetMaximum() << std::endl;
   pT_distr->SetParameter(0, 15);
    pT_distr->Draw();
}
