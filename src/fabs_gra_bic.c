
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

double myround(double a, int n)
{
    double s;
    if (a<0) {
        s=-a;
    }else{
        s=a;
    }
    s = s * pow(10.0, 1.0*n);
    s = s + 0.5;
    s = (int)s;
    s = s / pow(10.0, 1.0*n);
    if (a<0) {
        return -s;
    }else{
        return s;
    }
}

double L_h1(double u, double h) {
    if (u >= h) return u;
    else if (u <= -h) return -u;
    else return (3.0 * u * u) / (4.0 * h)
             - (u * u * u * u) / (8.0 * h * h * h)
             + (3.0 * h) / 8.0;
}

double L_h_prime(double u, double h) {
    if (u >= h) return 1.0;
    else if (u <= -h) return -1.0;
    else return (3.0 * u) / (2.0 * h)
             - (u * u * u) / (2.0 * h * h * h);
}

// difflaplace: Increment of Laplacian penalty
void loss(double *Loss, double *LLoss, double *y, double *xbeta, 
          double difflaplace, double *lambda2, double h, int *param)
{
    int n = param[0];
    int i, j;

    double total = 0.0;
    double diff_y, diff_xbeta, u;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            diff_y = y[i] - y[j];
            diff_xbeta = xbeta[i] - xbeta[j];
            u = diff_y - diff_xbeta;
            total += L_h1(u, h);
        }
    }

    *Loss = total / (n * (n - 1));

    *LLoss = *Loss + lambda2[0] * difflaplace;
}

// Laplacian matrix: L = D - A, i.e. Degree - Adjacency
// Dbeta: D * beta, Abeta: A * beta
void grad_with_Laplacian(double *der, double *y, double *x, double *xbeta,
    double *lambda2, double h, int *param, double *Dbeta, double *Abeta)
{
    int n = param[0];
    int p = param[1];
    int i, j, k;
    double diff_y, diff_xbeta, u, coeff;

    for (k = 0; k < p; k++) der[k] = 0.0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            diff_y = y[i] - y[j];
            diff_xbeta = xbeta[i] - xbeta[j];
            u = diff_y - diff_xbeta;
            coeff = -L_h_prime(u, h);  // 负号来自链式法则

            for (k = 0; k < p; k++) {
                der[k] += coeff * (x[i*p + k] - x[j*p + k]);
            }
        }
    }

    for (k = 0; k < p; k++) der[k] /= (n * (n - 1));

    if (lambda2[0] != 0) {
        for (i = 0; i < p; i++) der[i] += lambda2[0] * (Dbeta[i] - Abeta[i]);
    }
}

double calculate_EBIC(double *Loss, int *param, double gamma)
{
    int n = param[0];
    int p = param[1];
    int df = param[3];

    double ebic = 2.0 * n * (*Loss)
            + df * log((double)n)
            + 2.0 * gamma * lchoose((double)p, (double)df);

    return ebic;
}

// update sparse index
void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter)
{
    int i;
    s_i += *tt_a;
    s_j += *tt_a;
    for (i = 0; i < ns; ++i)
    {
        s_i[i] = act[i];
        s_j[i] = iter;
    }
    *tt_b = *tt_a;
    *tt_a += ns;
}

void update_Abeta(double *edgeweight, int *edgeweightLocRow, 
    int *edgeweightLocCol, double *Abeta, double value, int k, int *param)
{
    int i;

    for (i = 0; i < param[2]; ++i)
    {
        if ( edgeweightLocCol[i] == k ) {
            Abeta[ edgeweightLocRow[i] ] += edgeweight[i] * value;
        }
        if (edgeweightLocRow[i] == k) {
            Abeta[ edgeweightLocCol[i] ] += edgeweight[i] * value;
        }
    }
}

// take an initial step
void Initial(double *der, double *y, double *x, double *xbeta, int *param,
    double *lambda2, double h, double *Dbeta, double *Abeta, 
    double *weight, double eps, double *beta, double *Loss, double *LLoss, 
    double *degree, double *edge, int *edgeRow, int *edgeCol, double *lambda, 
    int *direction, int *active, double *bic, int *nactive, double *difflaplace)
{
    // der doesn't need initial value
    int i, k=0;
    int n = param[0], p = param[1];
    double temp, temp2, value;

    // calculate the derivative
    grad_with_Laplacian(der, y, x, xbeta, lambda2, h, param, Dbeta, Abeta);

    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i){
        value = fabs(der[i])/weight[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }

    // calculate increment, update xb, lambda, beta and loss_after
    value = eps/weight[k];
    if (der[k] > 0.0) value *= -1.0;
    *beta = value;
    
    // calculate Loss and LLoss for beta = 0.
    loss(&temp, &temp2, y, xbeta, *difflaplace, lambda2, h, param);

    // update xbeta, Loss, LLoss, difflaplace, Dbeta, Abeta, lambda
    for (i = 0; i < n; ++i) xbeta[i] = x[i*p + k] * value;
    Dbeta[k] = degree[k]*value;
    update_Abeta(edge, edgeRow, edgeCol, Abeta, value, k, param);
    *difflaplace = degree[k]*value*value/2.0;
    loss(Loss, LLoss, y, xbeta, *difflaplace, lambda2, h, param);
    *lambda = (temp2 - *LLoss)/eps;

    *direction = 1;
    param[3]   = 1;
    *active    = k;
    *bic       = calculate_EBIC(Loss, param, 1.0);
    *nactive   = 1;
}

int Backward(double *der, double *y, double *x, double *xbeta, int *param,
    double *lambda2, double h, double *Dbeta, double *Abeta, 
    double *weight, double eps, double *beta, double *Loss, double *LLoss, 
    double *degree, double *edge, int *edgeRow, int *edgeCol, 
    double *lambda, int *direction, int *active,double *bic, int *nactive, 
    double *beta0, double *xbtemp, double xi, double *difflaplace)
{
    int i, c, k;
    int n = param[0], p = param[1], ns = param[3];
    double temp, value;

    // update der
    grad_with_Laplacian(der, y, x, xbeta, lambda2, h, param, Dbeta, Abeta);

    // find a backward direction
    k = 0;
    c = active[0];
    temp = der[c] / weight[c];
    if (beta[0] < 0) temp *= -1.0;
    beta0[0]  = beta[0];
    for (i = 1; i < ns; ++i){
        c = active[i];
        value = der[c] / weight[c];
        if (beta[i] < 0) value *= -1.0;
        if (value > temp){
            temp = value;
            k = i;
        }
        beta0[i]  = beta[i];
    }

    // try a backward step.
    c = active[k];
    value = eps/weight[c];
    if (beta[k] > 0) value *= -1.0;
    beta0[k] += value;
    for (i = 0; i < n; ++i)
        xbtemp[i] = xbeta[i] + x[i*p+c]*value;
    temp = *difflaplace + degree[c]*value*value/2.0 + value*(Dbeta[c] - Abeta[c]);
    loss(Loss+1, LLoss+1, y, xbtemp, temp, lambda2, h, param);

    if (myround(LLoss[1] - LLoss[0] - (*lambda)*eps/weight[c], 6) < (-1.0*xi))
    {
        // adopt a backward step.
        // update direction, bic, nactive, lambda, difflaplace, Dbeta, Abeta, 
        // active, param[3], beta, xbeta, Loss, LLoss
        direction[1] = -1;
        bic[1]       = calculate_EBIC(Loss+1, param, 1);
        nactive[1]   = nactive[0];
        lambda[1]    = lambda[0];
        *difflaplace = temp;
        for (i = 0; i < n; ++i)
            xbeta[i] = xbtemp[i];
        update_Abeta(edge, edgeRow, edgeCol, Abeta, value, c, param);
        Dbeta[c] += degree[c]*value;
        // test if vanish
        if (fabs(beta0[k]) < 1.0*xi) {
            for (i = k; i < param[3]; ++i)
            {
                active[i] = active[i+1];
                beta0[i] = beta0[i+1];
            }
            param[3]--;
            nactive[1]--;
        }
        return 0;
    }
    beta0[k] -= value;
    return 1;
}

void Forward(double *der, double *y, double *x, double *xbeta, int *param,
    double *lambda2, double h, double *Dbeta, double *Abeta, double *weight, 
    double eps, double *beta, double *Loss, double *LLoss, double *degree, 
    double *edge, int *edgeRow, int *edgeCol, double *lambda, int *direction, 
    int *active, double *bic, int *nactive, double *difflaplace, double xi)
{
    // der doesn't need initial value
    int i, j, k=0;
    int n = param[0], p = param[1], ns = param[3];
    double temp, value;
    
    // d has calculated in Backward.
    // update der
    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i){
        value = fabs(der[i])/weight[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }

    // calculate increment
    value = eps/weight[k];
    if (der[k] > 0.0) value *= -1.0;
    // beta has been assigned in backward
    // update beta, active, nactive, param[3]
    if (k > active[ns-1])
    {
        active[ns] = k;
        beta[ns] = value;
        nactive[1] = nactive[0]+1;
        param[3]++;
    } else {
        for (i = 0; i < ns; ++i)
        {
            if (active[i] < k) continue;
            if (active[i] == k)
            {
                beta[i] += value;
                nactive[1] = nactive[0];
            } else {
                for (j = ns; j > i; --j)
                {
                    active[j] = active[j-1];
                    beta[j] = beta[j-1];
                }
                active[i] = k;
                beta[i] = value;
                param[3]++;
                nactive[1] = nactive[0]+1;
            }
            break;
        }
    }

    // update xbeta, Loss, LLoss, difflaplace, Dbeta, Abeta, lambda, direction, bic
    for (i = 0; i < n; ++i)
        xbeta[i] += x[i*p+k]*value;
    *difflaplace += degree[k]*value*value/2 + value*(Dbeta[k] - Abeta[k]);
    Dbeta[k] += degree[k]*value;
    update_Abeta(edge, edgeRow, edgeCol, Abeta, value, k, param);
    loss(Loss+1, LLoss+1, y, xbeta, difflaplace[0], lambda2, h, param);
    temp = (LLoss[0] - LLoss[1] - xi) * (weight[k] / eps);
    if (temp < *lambda) lambda[1] = temp;
    else lambda[1] = lambda[0];
    // test if vanish
    direction[1] = 1;
    bic[1]       = calculate_EBIC(Loss+1, param, 1);
}

void LFabs_single_lambda2(double *der, double *y, double *x, double *xbeta, int *param, 
    double *lambda2, double h, double *Dbeta, double *Abeta, double *weight, 
    double eps, double *beta, double *Loss, double *LLoss, double *degree, 
    double *edge, int *edgeRow, int *edgeCol, double *lambda, int *direction,
    int *active, double *bic, int*nactive, int iter, double *xbtemp,
    double xi, int max_s, int *sparse_i, int *sparse_j, double lam_m, int stoping)
{
    int i, k;
    double difflaplace = 0.0;
    
    // step 1: initial step (forward)
    Initial(der, y, x, xbeta, param, lambda2, h, Dbeta, Abeta,
            weight, eps, beta, Loss, LLoss, degree, edge, edgeRow, edgeCol,
            lambda, direction, active, bic, nactive, &difflaplace);
    int tt_act_b = 0;
    int tt_act_a = 0;
    usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], 0);
    
    // step 2: forward and backward
    for (i = 0; i < iter-1; ++i)
    {
        k = Backward(der, y, x, xbeta, param, lambda2, h, Dbeta, Abeta, 
                    weight, eps, beta+tt_act_b, Loss+i, LLoss+i, degree, 
                    edge, edgeRow, edgeCol, lambda+i, direction+i, active, 
                    bic+i, nactive+i, beta+tt_act_a, xbtemp, xi, &difflaplace);
        if (k) {
            Forward(der, y, x, xbeta, param, lambda2, h, Dbeta, Abeta,
                    weight, eps, beta+tt_act_a, Loss+i, LLoss+i, degree, 
                    edge, edgeRow, edgeCol, lambda+i, direction+i, active, 
                    bic+i, nactive+i, &difflaplace, xi);
        }
        usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], i+1);
        
        if ( stoping && (lambda[i+1] <= lambda[0] * lam_m) ) {
            iter = i+2;
            if (lambda[i+1] < 0)
            {
                iter--;
                tt_act_a -= param[3];
            }
            break;
        }
        if (param[3] > max_s) {
//            Rprintf("Warning! Max nonzero number is larger than predetermined threshold. Program ended early.\n");
            iter = i+2;
            break;
        }
        if (i == iter-2) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }
    
    param[4] = tt_act_a;
    param[5] = iter;
}

SEXP LFabs(SEXP Y, SEXP X, SEXP Weight, SEXP Epsilon, SEXP Lam_min,
           SEXP Xi, SEXP Stoping, SEXP Iter, SEXP Param, 
           SEXP Edge, SEXP EdgeRow, SEXP EdgeCol, SEXP Degree, 
           SEXP Max_S, SEXP Lambda2, SEXP H)
{
    int i, n, p, llambda2, nlambda2, iter, max_s;
    double *y, *x, *weight, h, eps, lam_m, xi, *degree, *edge, *lambda2;
    int *param, stoping, *edgeRow, *edgeCol;

    param     = INTEGER(Param);
    n         = param[0];
    p         = param[1];
    nlambda2  = param[6];
    y         = REAL(Y);
    x         = REAL(X);
    weight    = REAL(Weight);
    eps       = REAL(Epsilon)[0];
    lam_m     = REAL(Lam_min)[0];
    xi        = REAL(Xi)[0];
    stoping   = INTEGER(Stoping)[0];
    iter      = INTEGER(Iter)[0];
    max_s     = INTEGER(Max_S)[0];
    degree     = REAL(Degree);
    edge      = REAL(Edge);
    edgeRow   = INTEGER(EdgeRow);
    edgeCol   = INTEGER(EdgeCol);
    lambda2   = REAL(Lambda2);
    h         = REAL(H)[0];
    
    double *der, *beta, *lambda, *xbeta, *bic, *loss;
    double *lloss, *xbtemp, *Dbeta, *Abeta;
    int *direction, *nactive, *active;
    int *sparse_i, *sparse_j;

    beta      = (double*)malloc(sizeof(double)*iter*max_s);
    sparse_i  =    (int*)calloc(iter*max_s, sizeof(int));
    sparse_j  =    (int*)calloc(iter*max_s, sizeof(int));
    lambda    = (double*)calloc(iter, sizeof(double));
    direction =    (int*)malloc(sizeof(int)   *iter);
    bic       = (double*)malloc(sizeof(double)*iter);
    loss      = (double*)malloc(sizeof(double)*iter);
    lloss     = (double*)malloc(sizeof(double)*iter);
    nactive   =    (int*)malloc(sizeof(int)   *iter);
    xbtemp    = (double*)malloc(sizeof(double)  *n);
    xbeta     = (double*)calloc(n, sizeof(double)); // x^T beta.
    der       = (double*)calloc(p, sizeof(double)); // 1st order Taylor Formula.
    active    =    (int*)calloc(max_s+1, sizeof(int));
    Dbeta     = (double*)calloc(p, sizeof(double));
    Abeta     = (double*)calloc(p, sizeof(double));
    
    SEXP AllResult;
    PROTECT(AllResult = allocVector(VECSXP, nlambda2));
    
    for (llambda2 = 0; llambda2 < nlambda2; ++llambda2)
    {
        LFabs_single_lambda2(der, y, x, xbeta, param, lambda2+llambda2, h,
                         Dbeta, Abeta, weight, eps, beta, loss, lloss, degree, edge,
                         edgeRow, edgeCol, lambda, direction, active, bic, nactive, iter,
                         xbtemp, xi, max_s, sparse_i, sparse_j, lam_m, stoping);
        
        int tt_act_a = param[4];
        int niter = param[5];
        
        SEXP Beta, Lam, Drct, Loops, LLoss, BIC, Loss, Nactive, Idi, Idj, Result, R_names;
        char *names[10] = {"beta", "lambda", "direction", "iter", "bic", "loss",
            "losslaplace", "n_active", "index_i", "index_j"};
        PROTECT(Beta    = allocVector(REALSXP, tt_act_a));
        PROTECT(Idi     = allocVector(INTSXP,  tt_act_a));
        PROTECT(Idj     = allocVector(INTSXP,  tt_act_a));
        PROTECT(Lam     = allocVector(REALSXP, niter));
        PROTECT(BIC     = allocVector(REALSXP, niter));
        PROTECT(Loss    = allocVector(REALSXP, niter));
        PROTECT(LLoss   = allocVector(REALSXP, niter));
        PROTECT(Drct    = allocVector(INTSXP,  niter));
        PROTECT(Nactive = allocVector(INTSXP,  niter));
        PROTECT(Loops   = allocVector(INTSXP,  1));
        PROTECT(Result  = allocVector(VECSXP,  10));
        PROTECT(R_names = allocVector(STRSXP,  10));
        
        for(i = 0; i < 10; ++i) SET_STRING_ELT(R_names, i, mkChar(names[i]));
        INTEGER(Loops)[0] = niter;
        for (i = 0; i < tt_act_a; ++i)
        {
            REAL(Beta)[i]   = beta[i];
            INTEGER(Idi)[i] = sparse_i[i];
            INTEGER(Idj)[i] = sparse_j[i];
        }
        for (i = 0; i < niter; ++i)
        {
            REAL(BIC)[i]        = bic[i];
            REAL(Loss)[i]       = loss[i];
            REAL(LLoss)[i]      = lloss[i];
            REAL(Lam)[i]        = lambda[i];
            INTEGER(Drct)[i]    = direction[i];
            INTEGER(Nactive)[i] = nactive[i];
        }
        
        SET_VECTOR_ELT(Result, 0, Beta);
        SET_VECTOR_ELT(Result, 1, Lam);
        SET_VECTOR_ELT(Result, 2, Drct);
        SET_VECTOR_ELT(Result, 3, Loops);
        SET_VECTOR_ELT(Result, 4, BIC);
        SET_VECTOR_ELT(Result, 5, Loss);
        SET_VECTOR_ELT(Result, 6, LLoss);
        SET_VECTOR_ELT(Result, 7, Nactive);
        SET_VECTOR_ELT(Result, 8, Idi);
        SET_VECTOR_ELT(Result, 9, Idj);
        setAttrib(Result, R_NamesSymbol, R_names);
        
        SET_VECTOR_ELT(AllResult, llambda2, Result);
        
        if (llambda2 < nlambda2-1)
        {
            param[3] = 0;
            for (i = 0; i < tt_act_a; ++i)
            {
                beta[i]     = 0.0;
                sparse_i[i] = 0;
                sparse_j[i] = 0;
            }
            
            for (i = 0; i < niter; ++i)
            {
                bic[i]       = 0.0;
                loss[i]      = 0.0;
                lloss[i]     = 0.0;
                lambda[i]    = 0.0;
                direction[i] = 0;
                nactive[i]   = 0;
            }
            
            for (i = 0; i < n; ++i) xbeta[i] = 0.0; // Need to be initialized.
            for (i = 0; i < max_s+1; ++i) active[i] = 0;
            for (i = 0; i < p; ++i)
            {
                Dbeta[i] = 0.0;  // Need to be initialized.
                Abeta[i] = 0.0;
            }
        }
    }
    
    free(beta      );
    free(sparse_i  );
    free(sparse_j  );
    free(lambda    );
    free(direction );
    free(bic       );
    free(loss      );
    free(lloss     );
    free(nactive   );
    free(xbtemp    );
    free(xbeta     );
    free(der       );
    free(active    );
    free(Dbeta     );
    free(Abeta     );
    
    UNPROTECT(12*nlambda2+1);
    return AllResult;
}
