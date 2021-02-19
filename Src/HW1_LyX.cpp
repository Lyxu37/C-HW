#include "Homework1/Homework1.hpp"
#include "cfl/Interp.hpp"
#include <numeric>


using namespace prb;
using namespace std;
using namespace cfl;

cfl::Function volatilityCurve(double dLambda,double dInitialTime);

// Question 1

cfl::Function prb::forwardFX(double dSpotFX, const cfl::Function &rDomesticDiscount,
                             const cfl::Function &rForeignDiscount) {
    return dSpotFX*rForeignDiscount/rDomesticDiscount; // Function相除
    // ex first and save = save first and ex (covered interest rate parity)
    // 1 + domestic interest rate = 1/the spot rate *(1 + foreign interest rate)*forward rate
    // forward rate = the spot rate x (1 + domestic interest rate) / (1 + foreign interest rate)

}


//Question 3

cfl::Function prb::forwardLogLinearInterp(double dSpot, const std::vector<double> &rDeliveryTimes,
                                          const std::vector<double> &rForwardPrices, double dInitialTime) {
    // log can decrease oscillation and achieve good accuracy w/ low computational cost
    PRECONDITION(rDeliveryTimes.size() == rForwardPrices.size());
    PRECONDITION(!rDeliveryTimes.empty());
    PRECONDITION(rDeliveryTimes.front() > dInitialTime);
    PRECONDITION(std::equal(rDeliveryTimes.begin()+1,rDeliveryTimes.end(),rDeliveryTimes.begin(),
                            std::greater<double>())); //ensure t_1 < t_2 < t_3 < t_n...

    // initial time + forward time to be a new time vector
    std::vector<double> uTime(rDeliveryTimes);
    uTime.insert(uTime.begin(),dInitialTime);

    // logarithm of forward curves
    std::vector<double> LogFwd(uTime.size());
    LogFwd.front() = std::log(dSpot); // c.front() is equivalent to *c.begin()
    std::transform(rForwardPrices.begin(),rForwardPrices.end(),
                   LogFwd.begin()+1,[](double dX){return std::log(dX);});

    cfl::Interp uLinear = NInterp::linear();
    Function uLogFwdFunction = uLinear.interpolate(uTime.begin(),uTime.end(),LogFwd.begin());
    // Function u1LogFwdFunction = NInterp::linear().interpolate()

    return cfl::exp(uLogFwdFunction); //Function -- fwd curve
}


// Question 4

class VolatilityCurveHW: public cfl::IFunction
{
public:
    VolatilityCurveHW(double dLambda, double dInitialTime):
    m_dLambda(dLambda),m_dInitialTime(dInitialTime){}

    double operator()(double dT) const {
        PRECONDITION(belongs(dT));
        if (m_dLambda*(dT-m_dInitialTime) < cfl::c_dEps) {
            return 1.;
        }
        return std::sqrt((1.0-std::exp(-2.0*m_dLambda*(dT-m_dInitialTime)))
                /(2.0*m_dLambda*(dT-m_dInitialTime)));
    }

    bool belongs(double dT) const {
        return dT >= m_dInitialTime;
    }

private:
    double m_dLambda;
    double m_dInitialTime;
};

cfl::Function volatilityCurve(double dLambda,double dInitialTime){
    return cfl::Function (new VolatilityCurveHW(dLambda, dInitialTime));
}


cfl::Function prb::volatilityFitBlack(const std::vector<double> &rMaturities,
                                      const std::vector<double> &rVolatilities,
                                      double dLambda, double dInitialTime) {
    PRECONDITION(rMaturities.size() == rVolatilities.size());
    PRECONDITION(dLambda>0.);
    PRECONDITION(rMaturities.front() > dInitialTime);
    PRECONDITION(std::equal(rMaturities.begin()+1,rMaturities.end(),rMaturities.begin(),
                            std::greater<double>()));

    std::vector<double> model_sigma(rMaturities.size());
    std::transform(rMaturities.begin(),rMaturities.end(),model_sigma.begin(),
                   [dLambda,dInitialTime](double dX){
        return std::sqrt((1.0-std::exp(-2.0*dLambda*(dX-dInitialTime)))
                /(2.0*dLambda*(dX-dInitialTime)));});

    cfl::Function volatility_Curve = volatilityCurve(dLambda,dInitialTime);

    double dCov = std::inner_product(rVolatilities.begin(),rVolatilities.end(),
                                     model_sigma.begin(),0.0);
    double dVar = std::inner_product(rVolatilities.begin(),rVolatilities.end(),
                                     rVolatilities.begin(),0.0);
    double dKappa = dCov/dVar;
    return dKappa*volatility_Curve;
}


// Question 2

class ForwardCouponBond: public cfl::IFunction{
public:
    ForwardCouponBond(const cfl::Data::CashFlow &rBond, const cfl::Function &rDiscount,
                      double dInitialTime, bool bClean): m_ubond(rBond),m_uDiscount(rDiscount),
                      m_dInitialTime(dInitialTime),m_bClean(bClean){}

    double operator()(double dT) const {
        PRECONDITION(belongs(dT));
        int m = floor((dT - m_dInitialTime)/m_ubond.period) + 1;
        double sum_coupon(0.);
        for (unsigned int i(m);i<m_ubond.numberOfPayments;i++){
            sum_coupon += m_uDiscount(m_dInitialTime + m_ubond.period*i);
        }
        sum_coupon *= m_ubond.rate*m_ubond.period;
        double dirty = (sum_coupon + m_uDiscount(m_dInitialTime + m_ubond.numberOfPayments*m_ubond.period))
                /m_uDiscount(dT);
        if (m_bClean){
            double AccrTime = dT - (m_dInitialTime + m_ubond.period*(m-1));
            ASSERT(AccrTime>=0);
            ASSERT(AccrTime<m_ubond.period);
            dirty -= m_ubond.rate * AccrTime;
        }
        dirty *= m_ubond.notional;
        return dirty;
    }

    bool belongs(double dT) const {
        return ((dT >= m_dInitialTime) && (m_uDiscount.belongs(dT)) && (dT<=
        m_dInitialTime + m_ubond.numberOfPayments*m_ubond.period));
    }

private:
    cfl::Data::CashFlow m_ubond;
    cfl::Function m_uDiscount;
    double m_dInitialTime;
    bool m_bClean;
};

cfl::Function prb::forwardCouponBond(const cfl::Data::CashFlow &rBond, const cfl::Function &rDiscount,
                                     double dInitialTime, bool bClean) {
    return cfl::Function(new ForwardCouponBond(rBond,rDiscount,dInitialTime,bClean));
}