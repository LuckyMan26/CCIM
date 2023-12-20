﻿// Halide tutorial lesson 1: Getting started with Funcs, Vars, and Exprs

// This lesson demonstrates basic usage of Halide as a JIT compiler for imaging.

// On linux, you can compile and run it like so:
// g++ lesson_01*.cpp -g -I <path/to/Halide.h> -L <path/to/libHalide.so> -lHalide -lpthread -ldl -o lesson_01 -std=c++17
// LD_LIBRARY_PATH=<path/to/libHalide.so> ./lesson_01

// On os x:
// g++ lesson_01*.cpp -g -I <path/to/Halide.h> -L <path/to/libHalide.so> -lHalide -o lesson_01 -std=c++17
// DYLD_LIBRARY_PATH=<path/to/libHalide.dylib> ./lesson_01

// If you have the entire Halide source tree, you can also build it by
// running:
//    make tutorial_lesson_01_basics
// in a shell with the current directory at the top of the halide
// source tree.

// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"

// We'll also include stdio for printf.
#include <stdio.h>

Halide::Func matrix_power_p(Halide::Var x, Halide::Var y, Halide::Var d, Halide::Func& input, int exponent, int cols) {
   
    Halide::Func result;
    if (exponent == 0) {
        result(x, y,d) = select(x == y, 1, 0);
    }
    else if (exponent > 0) {
        int e = exponent - 1;
        result(x, y,d) = input((x - e) % cols, y,d);

    }

    else {
        int e = (exponent % cols + cols) - 1;
        result(x, y, d) = input((x + e) % cols, y,d);
    }


    return result;
}
Halide::Func matrix_multiplication(Halide::Var x, Halide::Var y, Halide::Var d,Halide::Func& L, Halide::Func& P) {
    Halide::Func res;
    Halide::RDom r(0, 4, "r");
    
    res(x, y,d) = Halide::sum(L(r, y,d) * P(x, r,d));
    std::cout << std::endl;
  
 
    return res;
}
Halide::Func matrix_power_q(Halide::Var x, Halide::Var y, Halide::Var d,Halide::Func& input, int exponent, int rows) {
   
    Halide::Func result;
    if (exponent == 0) {
        result(x, y,d) = select(
        (x==y),1,
            0
        );
    }
    else if (exponent > 0) {
        int e = exponent - 1;
        result(x, y,d) = input(x, (y + e) % rows,d);

    }

    else {
        int e = (exponent % rows + rows) - 1;
        result(x, y,d) = input(x, (y - e) % rows,d);

    }

  
  
    return result;
}

Halide::Func calculate_L(Halide::Func M, Halide::Var x, Halide::Var y, int m, int n, Halide::Var d, int s, int t) {
 
    Halide::Func L;
    L(x, y, d) = M(x, y, m - 1, n - 1, d) - M(x, y, m - 1, n - t - 1, d) - M(x, y, m - s - 1, n - 1, d) + M(x, y, m - s - 1, n - t - 1, d);

    return L;
}

Halide::Func calculate_G(Halide::Func M, Halide::Var x, Halide::Var y, int m, int n, Halide::Var d, int s, int t) {

    Halide::Func G;
    G(x, y, d) = M(x, y, m - 1, n - t - 1, d) - M(x, y, m - s - 1, n - 1, d);

    return G;
}

Halide::Func calculate_J(Halide::Func M, Halide::Var x, Halide::Var y, int m, int n, Halide::Var d, int s, int t) {
 
    Halide::Func J;
    J(x, y, d) = M(x, y, m - s - 1, n - t - 1, d);

    return J;
}

Halide::Func calculate_K(Halide::Func M, Halide::Var x, Halide::Var y, int m, int n, Halide::Var d, int s, int t) {
   
    Halide::Func K;
     K(x, y, d) = M(x, y, m - s- 1, n - 1, d) - M(x, y, m - s - 1, n - t- 1, d);

    return K;
}
int main(int argc, char** argv) {
  
    Halide::Var s, t, x, y, d;
    Halide::Func  X, B, M;

    //Here we set up matrix dimensions
    //It is size of target box
    int m, n;
    m = 4;
    n = 4;
    //it is size of full image
    int Mi, N;
    Mi = 4;
    N = 4;
   
 
    Halide::Func P;
    Halide::Func Q;

    P(x, y,d) = select(
        (x == 0 && y < m - 1), 0,
        (x == 0 && y == m - 1), 1,
        (y == m - 1 && x > 1), 0,
        (x > 0 && y < m - 1 && y == x - 1), 1,

        0
    );
    Q(x, y,d) = select(
        (y == 0 && x < n - 1), 0,
        (y == 0 && x == n - 1), 1,
        (y > 0 && x < n - 1 && y - 1 == x), 1,

        0
    );

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Halide::Func power_p = matrix_power_p(x,y,d,P, -i, Mi);
            Halide::Func power_q = matrix_power_q(x, y, d,Q, -j, N);
            B(x, y, s, t, d) = matrix_multiplication(x,y,d,power_p, power_q)(x, y, d);
        }
    }

   
   
   
  
    M(x, y, s, t, d) = 0;
    M(x, y, 0, 0, d) = B(x, y, 0, 0, d);


    for (int i = 1; i < m; i++) {
        M(x, y, i, 0, d) = B(x, y, i, 0, d) + M(x, y, i - 1, 0, d);
   }
    for (int i = 1; i < n; i++) {
        M(x, y, 0, i, d) = B(x, y, 0, i, d) + M(x, y, 0, i-1, d);
    }
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            M(x, y, i, j, d) = M(x, y, i - 1, j, d) + M(x, y, i, j - 1, d) - M(x, y, i - 1, j - 1, d) + B(x, y, i, j, d);
        }
    }
   
    
    Halide::Func L, G, K, J;
    //Here should be values of s,t(Z_s,t) 
    int S, T;

    S = 2;
    T = 2;
  
    L(x, y, d) = calculate_L(M,x,y,m,n,d,S,T)(x,y,d); 
    G(x, y, d) = calculate_G(M, x, y, m, n, d, S, T)(x, y, d);
    K(x, y, d) = calculate_K(M, x, y, m, n, d, S, T)(x, y, d);
    J(x, y, d) = calculate_J(M, x, y, m, n, d, S, T)(x, y, d);

  
    Halide::Func Z;
    Z(x,y,d) = L(x,y,d) + G(x, y, d) + K(x, y, d) + J(x, y, d);
    //Result is here
    Halide::Func res = Z.vectorize(d);

    
   
    return 0;
}
