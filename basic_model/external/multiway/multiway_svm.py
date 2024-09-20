import numpy as np
from sklearn.svm import LinearSVC

from ..multiway.multiway_model import MultiwayModel


class MultiwaySVM(MultiwayModel):
    def __init__(self, reg=1, penalty='l2', tolerance=1e-2, max_iter=100, J=252, K=799):
        super().__init__(reg, penalty, tolerance, max_iter, J, K)
        self.single_way_model_ = LinearSVC

    def set_params(self, **kwargs):
        super().set_params(**kwargs)

    def init_single_way_model(self, n):
        if n < self.J_ and n < self.K_:
            dual = True
        else:
            dual = False

        return self.single_way_model_(penalty=self.penalty, dual=dual, max_iter=10000, tol=1e-4)

    def compute_criterion(self, b0, b, X, y, reg):
        beta_X = b0 + X.dot(b)
        logv = np.mean(np.maximum(1 - y * beta_X, 0))
        if self.penalty == 'l2':
            logv = logv - reg * np.linalg.norm(b) ** 2 / 2
        elif self.penalty == 'l1':
            logv = logv - reg * np.sum(np.abs(b))
        return logv
