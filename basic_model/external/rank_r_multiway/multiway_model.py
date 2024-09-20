import numpy as np
import pandas as pd
from sklearn import metrics
import scipy as sc
from tensorly.decomposition import parafac


class MultiwayModel:
    def __init__(self, reg=1, penalty='l2', tolerance=1e-3, max_iter=100, J=252, K=799, rank=1, random_state=None,
                 pseudo_norm=False):
        self.reg_ = reg
        self.coef_ = None
        self.coef_K_ = None
        self.coef_J_ = None
        self.intercept_ = None
        self.penalty = penalty
        self.tolerance_ = tolerance
        self.max_iter_ = max_iter
        self.classes_ = None
        self.J_ = J
        self.K_ = K
        self.rank_ = rank
        self.logv_history_ = []
        self.single_way_model_ = None
        self.random_state_ = random_state
        self.pseudo_norm_ = pseudo_norm

    def set_params(self, **kwargs):
        self.reg_ = kwargs.get('reg', self.reg_)
        self.coef_ = kwargs.get('coef', self.coef_)
        self.penalty = kwargs.get('penalty', self.penalty)
        self.tolerance_ = kwargs.get('tolerance', self.tolerance_)
        self.max_iter_ = kwargs.get('max_iter', self.max_iter_)
        self.J_ = kwargs.get('J', self.J_)
        self.K_ = kwargs.get('K', self.K_)
        self.rank_ = kwargs.get('rank', self.rank_)
        self.random_state_ = kwargs.get('random_state', self.random_state_)
        self.pseudo_norm_ = kwargs.get('pseudo_norm', self.pseudo_norm_)

    def init_coef(self, direction, b=None, X=None, init='random'):
        if b is not None:
            return b
        if direction == 'J':
            N = self.J_
        else:
            N = self.K_

        if init == 'parafac' or init == 'random_parafac':
            _, factors = parafac(X, self.rank_)
            if direction == 'J':
                index = 1
            else:
                index = 2
            coef = factors[index]
            if init == 'random_parafac':
                coef += np.random.randn(N, self.rank_) * 1e-05

        elif init == 'mean_svd' or init == 'random_svd':
            X_mean = X.mean(axis=0)
            U, S, V_t = np.linalg.svd(X_mean)
            if direction == 'J':
                coef = U.T[:, :self.rank_]
            else:
                coef = V_t.T[:, :self.rank_]
            if init == 'random_svd':
                coef += np.random.randn(N, self.rank_) * 1e-05

        else:
            if self.random_state_ is not None:
                np.random.seed(self.random_state_)
            coef = np.random.randn(N, self.rank_)

        coef = coef / np.linalg.norm(coef, axis=0)
        return coef.reshape(N * self.rank_, 1, order='F')

    def init_single_way_model(self, n):
        return self.single_way_model_()

    def fit_single_way_model(self, model, X, y, nb_iter=10):
        if nb_iter < 1:
            model.set_params(tol=self.tolerance_ * 10)
        elif nb_iter < 2:
            model.set_params(tol=self.tolerance_)
        else:
            model.set_params(tol=self.tolerance_ * 0.1)
        model.fit(X, y)
        return model.coef_.copy().T, model.intercept_.copy()

    def compute_criterion(self, intercept, b, X, y, C):
        pass

    def get_inverse_of_regularisation_matrix(self, direction, b):
        if direction == 'J':
            N1, N2 = self.J_, self.K_
        else:
            N1, N2 = self.K_, self.J_
        R = np.zeros((self.rank_ * N1, self.rank_ * N1))

        if self.penalty == 'l2':
            if self.pseudo_norm_:
                for r in range(self.rank_):
                    R[r * N1:(r + 1) * N1, r * N1:(r + 1) * N1] = np.eye(N1) / (
                            np.linalg.norm(b[r * N2:(r + 1) * N2], ord=2) + 1e-09)
                return R
            B = b.reshape((N2, self.rank_), order='F')
            R = np.linalg.inv(sc.linalg.sqrtm(B.T.dot(B)).real)
            return np.kron(R, np.eye(N1))

        for r in range(self.rank_):
            R[r * N1:(r + 1) * N1, r * N1:(r + 1) * N1] = np.eye(N1) / (
                        np.linalg.norm(b[r * N2:(r + 1) * N2], ord=1) + 1e-09)
        return R

    def get_full_beta(self, direction, b1, b2):
        beta = np.zeros((self.J_ * self.K_, 1))
        if direction == 'J':
            b_J, b_K = b1, b2
        else:
            b_J, b_K = b2, b1
        for r in range(self.rank_):
            beta += np.kron(b_J[r * self.J_:(r + 1) * self.J_], b_K[r * self.K_:(r + 1) * self.K_])
        return beta

    def fit_step(self, direction, model, b, X, y, nb_iter=10):
        Z = []
        for r in range(self.rank_):
            if direction == 'J':
                Z.append(np.sum([b[k + r * self.K_] * X[:, :, k] for k in range(self.K_)], axis=0))
            else:
                Z.append(np.sum([b[j + r * self.J_] * X[:, j, :] for j in range(self.J_)], axis=0))
        Z = np.concatenate(Z, axis=1)
        R = self.get_inverse_of_regularisation_matrix(direction, b)
        Z = Z.dot(R)

        b2, intercept = self.fit_single_way_model(model, Z, y, nb_iter=nb_iter)
        self.logv_history_.append(self.compute_criterion(intercept.squeeze(), b2, Z, y, self.reg_))
        b2 = R.dot(b2)

        return b2, intercept

    def fit1D(self, X, y, b=None, init='random', first_direction='J'):
        n = len(y)
        if first_direction == 'J':
            second_direction = 'K'
        else:
            second_direction = 'J'
        b1 = None
        b2 = self.init_coef(second_direction, b=b, X=X, init=init)

        convergence = False
        nb_iter = 0
        intercept = np.array([0])

        model1 = self.init_single_way_model(n)
        model2 = self.init_single_way_model(n)

        while not convergence and nb_iter < self.max_iter_:
            b1, intercept = self.fit_step(first_direction, model1, b2, X, y, nb_iter)

            if np.linalg.norm(b1) < 1e-09:
                break

            b2, intercept = self.fit_step(second_direction, model2, b1, X, y, nb_iter)

            if np.linalg.norm(b2) < 1e-09:
                break

            nb_iter += 1
            if nb_iter >= self.max_iter_:
                print('Convergence warning, tolerance threshold has not been reached.')

            convergence = abs(self.logv_history_[-1] - self.logv_history_[-2]) / abs(
                self.logv_history_[-2]) < self.tolerance_

        B = self.get_full_beta(first_direction, b1, b2)
        b1, b2 = b1.reshape(-1, self.rank_, order='F'), b2.reshape(-1, self.rank_, order='F')
        if first_direction == 'J':
            return b2, b1, intercept.squeeze(), B
        return b1, b2, intercept.squeeze(), B

    def fit(self, X, y, b=None, init='random', first_direction='J'):
        if isinstance(X, pd.DataFrame):
            x = X.values.reshape((len(y), self.J_, self.K_))
        else:
            x = X.reshape((len(y), self.J_, self.K_))

        if isinstance(y, pd.Series):
            y = y.values

        classes = np.unique(y)
        if len(classes) > 2:
            self.classes_ = classes
            self.coef_K_, self.coef_J_, self.coef_, self.intercept_ = [], [], [], []

            for v in classes:
                y_v = (y == v).astype(int)
                b_K, b_J, b_0, B = self.fit1D(x, y_v, b=b, init=init, first_direction=first_direction)

                self.coef_K_.append(b_K)
                self.coef_J_.append(b_J)
                self.coef_.append(B)
                self.intercept_.append(b_0)

        else:
            self.coef_K_, self.coef_J_, self.intercept_, \
            self.coef_ = self.fit1D(x, y, b=b, init=init, first_direction=first_direction)

        return self

    def decision_function(self, X):
        if self.classes_ is None:
            if isinstance(X, pd.DataFrame):
                return X.values.dot(self.coef_) + self.intercept_
            return X.dot(self.coef_) + self.intercept_
        n = X.shape[0]
        c = len(self.classes_)
        scores = np.zeros((n, c))
        for i in range(c):
            if isinstance(X, pd.DataFrame):
                scores[:, i] = X.values.dot(self.coef_[i].squeeze()) + self.intercept_[i]
            scores[:, i] = X.dot(self.coef_[i].squeeze()) + self.intercept_[i]
        return scores

    def score(self, X, y):
        if self.classes_ is None:
            return metrics.roc_auc_score(y, self.decision_function(X))
        pred = np.argmax(self.decision_function(X), axis=1)
        y_pred = self.classes_[pred]
        return np.mean(y_pred == y)