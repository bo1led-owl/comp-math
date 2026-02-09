/*
machine epsilon:
minimal eps such that $1 + eps \ne 1$ and $1 + (eps / 2) = 1$

1. machine epsilon
2. fraction bits
3. 1,M*2^E, find E_max and E_min
4. compare 1, 1+(eps/2), 1+eps, 1+eps+(eps/2)
*/

#include <cassert>
#include <concepts>
#include <limits>
#include <print>

template <typename T>
void check(std::string_view name, T x, T expected) {
    std::println("{}: {} {} {}", name, x, (x == expected) ? "==" : "!=", expected);
}

template <typename T>
void compare(std::string_view xName, std::string_view yName, T x, T y) {
    if (x < y) {
        std::println("({}) < ({})", xName, yName);
    } else if (x > y) {
        std::println("({}) > ({})", xName, yName);
    } else {
        std::println("({}) = ({})", xName, yName);
    }
}

template <std::floating_point F>
void solve() {
    int fractionBits = 0;
    F eps = 1;
    while (F{1} + eps / 2 != F{1}) {
        eps /= 2;
        fractionBits++;
    }

    check("eps", eps, std::numeric_limits<F>::epsilon());
    check("fraction bits", fractionBits, std::numeric_limits<F>::digits - 1);

    int minExp = 0, maxExp = 0;
    {
        F v = 1;
        while (v != std::numeric_limits<F>::infinity()) {
            v *= 2;
            maxExp++;
        }
        maxExp--;

        v = 1;
        while (v != 0) {
            v /= 2;
            minExp--;
        }
        minExp++;
        minExp += fractionBits;
    }
    check("E_min", minExp, std::numeric_limits<F>::min_exponent - 1);
    check("E_max", maxExp, std::numeric_limits<F>::max_exponent - 1);

    std::array<F, 4> vals{1, 1 + eps / 2, 1 + eps, 1 + eps + eps / 2};
    std::array<std::string_view, 4> names{"1", "1 + e/2", "1 + e", "1 + e + e/2"};

    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i + 1; j < 4; ++j) {
            compare(names[i], names[j], vals[i], vals[j]);
        }
    }
}

int main() {
    std::println("float");
    solve<float>();
    std::println("\ndouble");
    solve<double>();
}
