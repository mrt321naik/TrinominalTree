import std;

// Here we build a double vector which will hold the tree.
std::vector<std::vector<double>> BOX(int n, int o)
{
  std::vector<std::vector<double>> res;
  std::vector<double> temp;

  for (int i {0}; i < n; ++i){
    temp.clear();
    for(int j {0}; j < o; ++j){
      temp.push_back(0.0);
    }
    res.push_back(temp);
  }
  return res;

}

// Print the output tree
void PRINTM(std::vector<std::vector<double>> x)
{
  for(auto & t : x)
  {
    for(auto & u : t)
    {
        std::print("{} \t", u);
    }
      std::cout << std::endl;
  }
}
// No price change factor
double one()
{
  return 1.0;
}
// Down price change factor
double down_factor (double up)
{
    return 1.0 / up; 
}

// No change Probability
double NoChange(double U, double D)
{
  return 1.0 - (U + D);
}
// Up price change factor
double up_factor(double v, double dt)
{
  return std::exp(v*std::sqrt(2.0*dt));
}

// Upper Probability
double up(double r, double q, double v, double dt)
{
  double a { std::exp((r-q)*dt/2.0)};
  double b { std::exp(-v*std::sqrt(dt/2.0))};
  double c { std::exp(v*std::sqrt(dt/2.0))};
  return std::pow((a - b)/(c - b), 2);
}

// Down Probability
double down(double r, double q, double v, double dt)
{
  double a { std::exp((r-q)*dt/2.0)};
  double b { std::exp(-v*std::sqrt(dt/2.0))};
  double c { std::exp(v*std::sqrt(dt/2.0))};
  return std::pow((c - a)/(c - b), 2);
}

// Calculates option price with trinominal tree
double PRICE(double S, double K, double r, double q, double v, double t, int nodes, std::string optype)
{ 
  // Calcualte measures
  int row = 4*nodes + 2;
  int col = nodes + 1;
  int split = row / 2 - 1;

  // Declare Tree
  std::vector<std::vector<double>> Tree = BOX(row, col);
  
  // Calculate delta time and factors
  double dt = t / (double) nodes;

  double UP = up_factor(v, dt);
  double DOWN = down_factor(UP);
  double M = one();

  double P_UP = up(r, q, v, dt);
  double P_DOWN = down(r, q, v, dt);
  double P_SAME = NoChange(P_UP,P_DOWN);
  

  // Set Initial stock price
  Tree[split][0] = S;
  
  // Forward propogation
  for(int j {};j < col; ++j ){

  
  for(int i {1}; i < col - j; ++i)
  {
    Tree[split - 2*i][j + i] = Tree[split - 2*(i-1)][j + i - 1]*UP;
    Tree[split + 2*i][j + i] = Tree[split + 2*(i-1)][j + i - 1]*DOWN;
    Tree[split][j + i] = Tree[split][j + i - 1]*M;
  }
}

  // Compute payoffs
  for(std::size_t i {1}; i < row; ++i)
  {
    if(i % 2 !=0){
      if(optype == "call"){
      Tree[i][col-1] = std::fmax(Tree[i - 1][col - 1] - K, 0.0);
    } else {
        Tree [i][col -1] = std::fmax(K - Tree[i - 1][col - 1], 0.0);
      }
    }
  }

  // Discount to compute option price

   double E = std::exp(-r*dt); // Discount factor 
   double A, B, C; 

  int f {3};
  for (int i {col - 2}; i>=0; --i){
    for(int j {f}; j <= row - f; ++j){
      if (j % 2 != 0){
        A = Tree[j - 2][i + 1];
        B = Tree[j][i + 1];
        C = Tree[j + 2][i + 1];
        if (optype == "call"){
      Tree[j][i] = std::fmax(E*(A*P_UP + B*P_SAME + C*P_DOWN), Tree[j - 1][i] - K ); 
        } else {
          Tree[j][i] = std::fmax(E*(A*P_UP + B*P_SAME + C*P_DOWN), K - Tree[j -1][i]);
        }
      }
    }

    f += 2;
  }
  

  //PRINTM(Tree);

  return Tree[split + 1][0]; 
  
}


int main()
{     
  double S {100};     // asset price
  double K {105};      // Strike price
  double r {0.05};    // Interest rate 
  double q {0.01};    // Divident yeild 
  double v {0.2};     // Volitality 
  double t {30.0/365.0};    // experation date
  int nodes {16};            // number of nodes
  std::string optype {"call"}; // Option type here is put

    //PRICE(S, K, r, q, v, t, nodes, optype);

  // IMPLIED VOLITALITY CALCULATOIN 
  double diff = 0, vega = 0, option_price = 2.50;
  double v0 = 0.1, v1 = 0.99, dv = 0.01;

  while(true){
    diff = PRICE(S, K, r, q, v0, t, nodes, optype) - option_price;
    vega = (PRICE(S, K, r, q, v0+dv, t, nodes, optype) - PRICE(S, K, r, q, v0-dv, t, nodes, optype))/(2.0*dv);

    v1 = v0 - diff / vega;
    if(std::abs(v1 - v0) < 0.0001) {
      break;
    }

    v0 = v1;

  }
    std::println("Option price: {}", option_price);
    std::println("Implied Volitality: {}", v1);
    std::println("Check Price {}",PRICE(S, K, r, q, v1, t, nodes, optype));

    // Calculate Greeks
    double delta, gamma, theta, rho;
    double dS = 0.01*S, dT = 1.0/365.0, dV = 0.01, dR = 0.01;

    delta = (PRICE(S+dS, K, r, q, v, t, nodes, optype) - PRICE(S-dS, K, r, q, v, t, nodes, optype)) / (2.0*dS);
    gamma = (PRICE(S+dS, K, r, q, v, t, nodes, optype) - 2.0*PRICE(S, K, r, q, v, t, nodes, optype) + PRICE(S-dS, K, r, q, v, t, nodes, optype)) / std::pow(dS, 2);
    theta = -(PRICE(S, K, r, q, v, t+dT, nodes, optype) - PRICE(S, K, r, q, v, t, nodes, optype)) / dT;
    vega = (PRICE(S, K, r, q, v+dV, t, nodes, optype) - PRICE(S, K, r, q, v-dV, t, nodes, optype))/(2.0*dv);
    rho = (PRICE(S, K, r+dR, q, v, t, nodes, optype) - PRICE(S, K, r-dR, q, v, t, nodes, optype)) / (2.0*dR);
    
    vega /= 100;
    rho  /= 100;

  std::println("");
  std::println("Delta: {}", delta);
  std::println("Gamma: {}", gamma);
  std::println("Theta: {}", theta);
  std::println("Vega: {}", vega);
  std::println("Rho: {}", rho);
  std::println("");

  return 0;
}
