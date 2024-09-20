import numpy as np
import pandas as pd
from sklearn import metrics


class MultiwayModel:
    def __init__(self, reg=1, penalty='l2', tolerance=1e-3, max_iter=100, J=252, K=799, random_state=None):
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
        self.logv_history_ = []
        self.single_way_model_ = None
        self.random_state_ = random_state

    def set_params(self, **kwargs):
        self.reg_ = kwargs.get('reg', self.reg_)
        self.coef_ = kwargs.get('coef', self.coef_)
        self.penalty = kwargs.get('penalty', self.penalty)
        self.tolerance_ = kwargs.get('tolerance', self.tolerance_)
        self.max_iter_ = kwargs.get('max_iter', self.max_iter_)
        self.J_ = kwargs.get('J', self.J_)
        self.K_ = kwargs.get('K', self.K_)
        self.random_state_ = kwargs.get('random_state', self.random_state_)

    def init_coef(self, direction, b=None, X=None, init='random'):
        if b is not None:
            return b
        if direction == 'J':
            N = self.J_
        else:
            N = self.K_
        if init == 'random':
            if self.random_state_ is not None:
                np.random.seed(self.random_state_)
            coef = np.random.randn(N, 1)
            return coef / np.linalg.norm(coef)

        if init == 'parietal_mean':
            if direction == 'J':
                coef = X.mean(axis=(0, 2)).reshape(N, 1)
            else:
                coef = X.mean(axis=(0, 1)).reshape(N, 1)
            return coef / np.linalg.norm(coef)

        if init == 'normalized_constant':
            return np.ones((N, 1)) / N

    def init_single_way_model(self, n):
        return self.single_way_model_()

    def fit_single_way_model(self, model, intercept, b, C, X, y, nb_iter=10):
        if b is not None:
            model.coef_ = b.copy().T
            model.intercept_ = intercept.copy()

        model.set_params(C=C)
        if nb_iter < 1:
            model.set_params(tol=self.tolerance_ * 10)
        elif nb_iter < 2:
            model.set_params(tol=self.tolerance_)
        else:
            model.set_params(tol=self.tolerance_ * 0.1)
        model.fit(X, y)
        return model.coef_.copy().T, model.intercept_.copy()

    def get_inverse_of_regularisation_strength(self, b):
        if self.penalty == 'l2':
            return 1 / (self.reg_ * np.linalg.norm(b) ** 2)
        return 1 / (self.reg_ * np.sum(np.abs(b)))

    def compute_criterion(self, intercept, b, X, y, C):
        pass

    def fit_step(self, direction, model, intercept, b1, b2, X, y, nb_iter=10):
        if direction == 'J':
            Z = np.sum([b1[k] * X[:, :, k] for k in range(self.K_)], axis=0)
        else:
            Z = np.sum([b1[j] * X[:, j, :] for j in range(self.J_)], axis=0)

        C = self.get_inverse_of_regularisation_strength(b1)
        b2, intercept = self.fit_single_way_model(model, intercept, b2, C, Z, y, nb_iter=nb_iter)
        self.logv_history_.append(self.compute_criterion(intercept.squeeze(), b2, Z, y, 1 / C))

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
            b1, intercept = self.fit_step(first_direction, model1, intercept, b2, b1, X, y, nb_iter)

            if np.linalg.norm(b1) < 1e-09:
                break

            b2, intercept = self.fit_step(second_direction, model2, intercept, b1, b2, X, y, nb_iter)

            if np.linalg.norm(b2) < 1e-09:
                break

            nb_iter += 1
            if nb_iter >= self.max_iter_:
                print('Convergence warning, tolerance threshold has not been reached.')

            convergence = abs(self.logv_history_[-1] - self.logv_history_[-2]) / abs(
                self.logv_history_[-2]) < self.tolerance_

        b1, b2 = b1.squeeze(), b2.squeeze()
        if first_direction == 'J':
            return b2 / np.linalg.norm(b2), b1 * np.linalg.norm(b2), intercept.squeeze(), np.kron(b1, b2)
        return b1 / np.linalg.norm(b1), b2 * np.linalg.norm(b1), intercept.squeeze(), np.kron(b2, b1)

    def fit(self, X, y, b=None, init='random', first_direction='J'):
        if isinstance(X, pd.DataFrame):
            x = X.values.reshape(len(y), self.J_, self.K_)
        else:
            x = X.reshape(len(y), self.J_, self.K_)

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
