import numpy as np
from sklearn.linear_model import LogisticRegression

from ..rank_r_multiway.multiway_model import MultiwayModel


class MultiwayLogisticRegression(MultiwayModel):
    def __init__(self, reg=1, penalty='l2', tolerance=1e-3, max_iter=100, J=252, K=799, solver='saga', n_jobs=1,
                 rank=1, random_state=None, pseudo_norm=False):
        super().__init__(reg, penalty, tolerance, max_iter, J, K, rank, random_state, pseudo_norm)
        self.solver_ = solver
        self.n_jobs_ = n_jobs
        self.single_way_model_ = LogisticRegression

    def set_params(self, **kwargs):
        super().set_params(**kwargs)
        self.n_jobs_ = kwargs.get('n_jobs', self.n_jobs_)
        self.solver_ = kwargs.get('solver', self.solver_)

    def init_single_way_model(self, n):
        if self.solver_ == 'liblinear' and n < self.J_ and n < self.K_:
            dual = True
        else:
            dual = False

        return self.single_way_model_(solver=self.solver_, penalty=self.penalty, dual=dual, C=1/self.reg_,
                                      max_iter=10000, n_jobs=self.n_jobs_, tol=self.tolerance_, warm_start=True,
                                      random_state=self.random_state_)

    def compute_criterion(self, b0, b, X, y, reg):
        beta_X = b0 + X.dot(b)
        logv = y.reshape(1, -1).dot(beta_X).squeeze() - np.sum(np.log(1 + np.exp(beta_X)))
        if self.penalty == 'l2':
            logv = logv - reg * np.linalg.norm(b) ** 2 / 2
        elif self.penalty == 'l1':
            logv = logv - reg * np.sum(np.abs(b))
        return logv
