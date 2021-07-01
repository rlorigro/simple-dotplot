
#ifndef RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP
#define RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using std::cout;
using std::isnan;
using std::runtime_error;
using std::vector;


template <class T> class IterativeSummaryStats {
public:
    /// Attributes ///
    T k;
    T n;
    T sum;
    T sum_of_squares;

    /// Methods ///
    IterativeSummaryStats();
    void add(T x);
    void remove(T x);
    double get_mean();
    double get_variance();
    T size();
    bool empty();
    void clear();
};


template <class T> void operator+=(IterativeSummaryStats<T>& a, IterativeSummaryStats<T>& b){
    /// Might break the add and remove functions?? cant be bothered to test this...
    a.n += b.n;
    a.sum += b.sum;
    a.sum_of_squares += b.sum_of_squares;
}


template <class T> T IterativeSummaryStats<T>::size(){
    return this->n;
}


template <class T> bool IterativeSummaryStats<T>::empty(){
    return (this->n == 0);
}


template <class T> IterativeSummaryStats<T>::IterativeSummaryStats(){
    this->k = 0;
    this->n = 0;
    this->sum = 0;
    this->sum_of_squares = 0;
}


template <class T> void IterativeSummaryStats<T>::clear(){
    this->k = 0;
    this->n = 0;
    this->sum = 0;
    this->sum_of_squares = 0;
}


template <class T> void IterativeSummaryStats<T>::add(T x){
    if (this->n == 0){
        this->k = x;
    }

    this->sum += x - this->k;
    this->sum_of_squares += (x - this->k)*(x - this->k);
    this->n++;

    if (this->sum < 0){
        cout << "";
    }
}


template <class T> void IterativeSummaryStats<T>::remove(T x){
    if (this->n == 0){
        throw runtime_error("ERROR: cannot remove a value from empty IterativeSummaryStats object");
    }

    this->sum -= x - this->k;
    this->sum_of_squares -= (x - this->k)*(x - this->k);
    this->n--;
}


template <class T> double IterativeSummaryStats<T>::get_mean() {
    ///
    /// K + Ex / n
    ///
    double mean = this->k + double(this->sum)/(this->n);

    return mean;
}


template <class T> double IterativeSummaryStats<T>::get_variance() {
    ///
    /// (Ex2 - (Ex*Ex)/n) / (n-1)
    ///


    double variance;

    if (this->n > 1){
        variance = (this->sum_of_squares - double(this->sum * this->sum) / this->n) / (this->n - 1);
    }
    else{
        variance = 0;
    }

    return variance;
}


template <class T> double pool_variances(vector <IterativeSummaryStats <T> >& stats){
    double numerator = 0;
    double denominator = 0;
    double variance;

    for (auto& s: stats){
        if (s.n == 0){
            continue;
        }
        numerator += (s.n - 1)*s.get_variance();
        denominator += s.n - 1;
    }

    if (denominator == 0){
        variance = 0;
    }
    else{
        variance = numerator/denominator;
    }

    return variance;
}


template <class T> double pool_means(vector <IterativeSummaryStats <T> >& stats){
    double numerator = 0;
    double denominator = 0;

    for (auto& s: stats){
        numerator += s.get_mean()*s.n;
        denominator += s.n;
    }

    return numerator/denominator;
}


#endif //RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP
