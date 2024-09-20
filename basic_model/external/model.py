import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
from time import time

from sklearn.model_selection import train_test_split, ParameterGrid, StratifiedKFold
from sklearn import metrics


class DataModel:

    def __init__(self, normalized=True, reduced=False, base_path='../../data/Python/', time_range=None, elec_range=None,
                 ind_range=None, pho_range=None):
        if normalized:
            if reduced:
                file_name = 'pho_normalized_reduced_array.npy'
            else:
                file_name = 'pho_normalized_array.npy'
        else:
            if reduced:
                file_name = 'pho_reduced_array.npy'
            else:
                file_name = 'pho_not_normalized_array.npy'
        np.load(base_path + file_name)
        self.file_path_ = base_path + file_name

        if time_range is None:
            if reduced:
                time_range = np.linspace(0, 998, 500)
            else:
                time_range = np.linspace(-200 + normalized * 2, 1398, 800 - normalized)

        if elec_range is None:
            elec_range = [f'E{i}' for i in range(1, 125)] + [f'E{i}' for i in range(129, 257)]

        if ind_range is None:
            ind_range = [f'evoked_df_S{"{:02d}".format((i))}-st-epo' for i in range(4, 42) if
                         i not in [6, 10, 11, 13, 16, 19, 20, 21, 25, 28, 30, 33, 36]]

        if pho_range is None:
            pho_range = ['bi_f', 'bi_m', 'bo_f', 'bo_m', 'di_f', 'di_m', 'do_f', 'do_m', 'gi_m', 'gi_f', 'go_m', 'go_f',
                         'ii_f', 'ii_m',
                         'mi_f', 'mi_m', 'mo_f', 'mo_m', 'ni_f', 'ni_m', 'no_m', 'no_f', 'gni_m', 'gni_f', 'gno_m',
                         'gno_f', 'oo_m', 'oo_f']

        self.time_range_ = time_range
        self.elec_range_ = elec_range
        self.ind_range_ = ind_range
        self.pho_range_ = pho_range

        self.multi_col_index_ = pd.MultiIndex.from_product([elec_range, time_range], names=['E', 'time'])
        self.multi_row_index_ = pd.MultiIndex.from_product([pho_range, ind_range], names=['pho', 'ind'])

        self.bilabials_ = ['bi_f', 'bi_m', 'bo_f', 'bo_m', 'mi_f', 'mi_m', 'mo_f', 'mo_m']
        self.alveolars_ = ['di_f', 'di_m', 'do_f', 'do_m', 'ni_f', 'ni_m', 'no_m', 'no_f']
        self.velars_ = ['gi_m', 'gi_f', 'go_m', 'go_f', 'gni_m', 'gni_f', 'gno_m', 'gno_f']

        self.obstruents_ = ['bi_f', 'bi_m', 'bo_f', 'bo_m', 'di_f', 'di_m', 'do_f', 'do_m', 'gi_m', 'gi_f', 'go_m',
                            'go_f']
        self.sonorants_ = ['mi_f', 'mi_m', 'mo_f', 'mo_m', 'ni_f', 'ni_m', 'no_m', 'no_f', 'gni_m', 'gni_f', 'gno_m',
                           'gno_f']
        self.normalized_ = normalized
        self.reduced_ = reduced

    def get_df(self):
        data = np.load(self.file_path_)
        return pd.DataFrame(
            data=data.reshape(-1, len(self.ind_range_) * len(self.pho_range_)).T,
            index=self.multi_row_index_,
            columns=self.multi_col_index_
        )

    def label_data(self, problem='manner', df=None):
        def affect_value(row, arrays):
            for i, array in enumerate(arrays):
                if row.name[0] in array:
                    return i
            return -1

        arrays = None
        if problem == 'manner':
            arrays = [self.obstruents_, self.sonorants_]
        if problem == 'place':
            arrays = [self.bilabials_, self.alveolars_, self.velars_]

        if df is None:
            return self.get_df().apply(affect_value, 1, args=[arrays])
        return df.apply(affect_value, 1, args=[arrays])

    def get_gender_and_vowel_features(self, df=None):
        if df is None:
            df = self.get_df()
        vowel = df.apply(lambda x: int(x.name[0][-3] == 'o'), 1)
        gender = df.apply(lambda x: int(x.name[0][-1] == 'm'), 1)
        return vowel, gender

    def train_test_split(self, test_size=0.2, problem='manner', remove=True, random_state=0):
        y = self.label_data(problem)
        if remove:
            y = y[y >= 0]

        X = self.get_df().reindex(y.index)
        ind_range_train, ind_range_test = train_test_split(self.ind_range_, test_size=test_size,
                                                           random_state=random_state)
        X_train, X_test = X.loc[pd.IndexSlice[:, ind_range_train], :], X.loc[pd.IndexSlice[:, ind_range_test], :]
        y_train, y_test = y.loc[:, ind_range_train], y.loc[:, ind_range_test]

        return X_train, X_test, y_train, y_test

    def CV_score(self, X, y, model, random_state=None, n_folds=5, problem='manner'):
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_state)
        val_scores = []
        train_scores = []
        for train_index, test_index in skf.split(X, y):
            X_cv_train, X_cv_test = X.reindex(X.index[train_index]), X.reindex(X.index[test_index])
            y_cv_train, y_cv_test = y.reindex(y.index[train_index]), y.reindex(y.index[test_index])

            model.fit(X_cv_train, y_cv_train)

            if problem == 'manner':
                if hasattr(model, 'decision_function'):
                    y_val_scores = model.decision_function(X_cv_test)
                    val_score = metrics.roc_auc_score(y_cv_test, y_val_scores)
                    val_scores.append(val_score)

                    y_train_scores = model.decision_function(X_cv_train)
                    train_score = metrics.roc_auc_score(y_cv_train, y_train_scores)
                    train_scores.append(train_score)
                else:
                    val_score = model.score(X_cv_test, y_cv_test)
                    val_scores.append(val_score)

                    train_score = model.score(X_cv_train, y_cv_train)
                    train_scores.append(train_score)

            if problem == 'place':
                val_score = model.score(X_cv_test, y_cv_test)
                val_scores.append(val_score)

                train_score = model.score(X_cv_train, y_cv_train)
                train_scores.append(train_score)

        return sum(val_scores) / len(val_scores), sum(train_scores) / len(train_scores)

    def grid_search_CV(self, model, param_grid, n_folds=5, problem='manner'):
        X_train, _, y_train, _ = self.train_test_split(problem=problem)
        results = {}
        for param_set in list(ParameterGrid(param_grid)):
            model.set_params(**param_set)
            val_score, train_score = self.CV_score(X_train, y_train, model, n_folds=n_folds, problem=problem)
            results[tuple(param_set.items())] = (val_score, train_score)
            print('param_set:', tuple(param_set.items()), 'validation score:', val_score, 'train score:', train_score)

        best_param_set = max(results, key=results.get)

        return best_param_set, results

    def score_leave_one_infant_out(self, model, problem='manner', verbose=True, init='normalized_constant',
                                   first_direction='J'):
        test_scores = []
        train_scores = []

        y = self.label_data(problem)
        y = y[y >= 0]

        for ind in self.ind_range_:
            df = self.get_df()
            df = df.reindex(y.index)

            X_train = df.drop(ind, level=1, axis=0)
            y_train = y.drop(ind, level=1, axis=0)
            X_test = df.xs(ind, level=1, drop_level=False)
            y_test = y.xs(ind, level=1, drop_level=False)

            df = None

            if 'Multiway' in type(model).__name__:
                model.fit(X_train, y_train, init=init, first_direction=first_direction)
            else:
                model.fit(X_train, y_train)
            if problem == 'manner':
                if hasattr(model, 'decision_function'):
                    y_test_scores = model.decision_function(X_test)
                    test_score = metrics.roc_auc_score(y_test, y_test_scores)
                    test_scores.append(test_score)

                    y_train_scores = model.decision_function(X_train)
                    train_score = metrics.roc_auc_score(y_train, y_train_scores)
                    train_scores.append(train_score)
                else:
                    test_score = model.score(X_test, y_test)
                    test_scores.append(test_score)

                    train_score = model.score(X_train, y_train)
                    train_scores.append(train_score)

            if problem == 'place':
                test_score = model.score(X_test, y_test)
                test_scores.append(test_score)

                train_score = model.score(X_train, y_train)
                train_scores.append(train_score)

            if verbose:
                print(ind, 'done \n', test_score, train_score)

        # return sum(test_scores) / len(test_scores), sum(train_scores) / len(train_scores)
        return test_scores, train_scores

    def grid_search_leave_one_infant_out(self, m, param_grid, problem='manner', write_to_file=True, verbose=2,
                                         init='normalized_constant', first_direction='J', base_path=''):
        results = {}
        for param_set in list(ParameterGrid(param_grid)):
            model = deepcopy(m)
            model.set_params(**param_set)
            start_time = time()
            val_scores, train_scores = self.score_leave_one_infant_out(model, problem=problem, verbose=verbose > 1,
                                                                       init=init, first_direction=first_direction)
            elapsed_time = time() - start_time
            val_mean, train_mean, val_var, train_var = np.mean(val_scores), np.mean(train_scores), np.var(
                val_scores), np.var(train_scores)
            results[tuple(param_set.items())] = (val_mean, train_mean)
            if verbose > 0:
                print('param_set:', tuple(param_set.items()), 'validation score:', val_mean, 'train score:', train_mean)
            if write_to_file:
                file_name = f'{base_path}{model.__class__.__name__}_{problem}_{model.penalty}_loio_test_results.csv'
                with open(file_name, 'a') as file:
                    file.write(','.join([str(val_mean), str(val_var), *[str(s) for s in val_scores],
                                         *[str(p) for p in param_set.values()], str(elapsed_time)]))
                    file.write('\n')
                file_name = f'{base_path}{model.__class__.__name__}_{problem}_{model.penalty}_loio_train_results.csv'
                with open(file_name, 'a') as file:
                    file.write(','.join([str(train_mean), str(train_var), *[str(s) for s in train_scores],
                                         *[str(p) for p in param_set.values()], str(elapsed_time)]))
                    file.write('\n')

        best_param_set = max(results, key=results.get)

        return best_param_set, results

    def eval_score(self, model, problem='manner'):
        X_train, X_test, y_train, y_test = self.train_test_split(problem=problem)
        model.fit(X_train, y_train)

        score = 0
        if problem == 'manner':
            if hasattr(model, 'decision_function'):
                y_scores = model.decision_function(X_test)
                fpr, tpr, _ = metrics.roc_curve(y_test, y_scores)
                score = metrics.roc_auc_score(y_test, y_scores)

                plt.plot(fpr, tpr)
                plt.show()
                print(f'AUC score: {score}')
            else:
                score = model.score(X_test, y_test)
                print(f'{score * 100} % of correct classifications')

        if problem == 'place':
            score = model.score(X_test, y_test)
            print(f'{score * 100} % of correct classifications')

        return model, score, X_train, y_train, X_test, y_test

    def compute_bootstrap_weights(self, model, problem='manner', B=300, B_size=25, write_to_File=True, init='random',
                                  first_direction='J', base_path=''):
        coefs, coefs_J, coefs_K, intercepts = [], [], [], []
        for _ in range(B):
            indice = np.random.choice(self.ind_range_, size=B_size, replace=True)
            X = self.get_df().reindex(pd.MultiIndex.from_product([self.multi_row_index_.levels[0], indice]))
            y = self.label_data(problem=problem, df=X)
            X = X[y >= 0]
            y = y[y >= 0]
            model_b = deepcopy(model)
            if 'Multiway' in model.__class__.__name__:
                model_b.fit(X, y, init=init, first_direction=first_direction)
            else:
                model_b.fit(X, y)
            if write_to_File:
                time_code = str(int(time()))
                if 'Multiway' in model.__class__.__name__:
                    np.save(f'{base_path}{time_code}_{model.__class__.__name__}_{model.penalty}_J_weights.npy',
                            model_b.coef_J_)
                    np.save(f'{base_path}{time_code}_{model.__class__.__name__}_{model.penalty}_K_weights.npy',
                            model_b.coef_K_)
                else:
                    np.save(f'{base_path}{time_code}_{model.__class__.__name__}_{model.penalty}_weights.npy',
                            model_b.coef_)
                np.save(f'{base_path}{time_code}_{model.__class__.__name__}_{model.penalty}_intercept.npy',
                        model_b.intercept_)
            if 'Multiway' in model.__class__.__name__:
                coefs_J.append(model_b.coef_J_)
                coefs_K.append(model_b.coef_K_)
            else:
                coefs.append(model.coef_)
            intercepts.append(model.intercept_)
        if 'Multiway' in model.__class__.__name__:
            return coefs_J, coefs_K, intercepts
        else:
            return coefs, intercepts

    def evaluate_weight_stability(self, model, multiway, problem='manner', B=300, B_size=25, sparse=False,
                                  write_to_file=True, init='random', first_direction='J'):
        def fleiss_kappa(coefs, B):
            coef_matrix = np.array([B - np.sum(np.array(coefs) != 0, axis=0), np.sum(np.array(coefs) != 0, axis=0)])
            P_i = 1 / (B * (B - 1)) * (coef_matrix[0] * (coef_matrix[0] - 1) + coef_matrix[1] * (coef_matrix[1] - 1))
            P_bar = np.mean(P_i)
            P_j = 1 / (len(coefs[0]) * B) * np.sum(coef_matrix, axis=1)
            P_e = np.sum(P_j ** 2)
            kappa = (P_bar - P_e) / (1 - P_e)
            freq = coef_matrix[1] / B

            return P_i, kappa, freq

        if problem == 'place':
            if multiway:
                coefs_J_bilabials, coefs_J_alveolars, coefs_J_velars = [], [], []
                coefs_K_bilabials, coefs_K_alveolars, coefs_K_velars = [], [], []
            else:
                coefs_bilabials, coefs_alveolars, coefs_velars = [], [], []
        else:
            if multiway:
                coefs_J, coefs_K = [], []
            else:
                coefs = []
        for _ in range(B):
            indice = np.random.choice(self.ind_range_, size=B_size, replace=True)
            X = self.get_df().reindex(pd.MultiIndex.from_product([self.multi_row_index_.levels[0], indice]))
            y = self.label_data(problem=problem, df=X)
            X = X[y >= 0]
            y = y[y >= 0]

            model_b = deepcopy(model)
            if multiway:
                model_b.fit(X, y, init=init, first_direction=first_direction)
            else:
                model_b.fit(X, y)

            if problem == 'place':
                if multiway:
                    coefs_J_bilabials.append(model_b.coef_J_[0])
                    coefs_J_alveolars.append(model_b.coef_J_[1])
                    coefs_J_velars.append(model_b.coef_J_[2])
                    coefs_K_bilabials.append(model_b.coef_K_[0])
                    coefs_K_alveolars.append(model_b.coef_K_[1])
                    coefs_K_velars.append(model_b.coef_K_[2])
                else:
                    coefs_bilabials.append(model_b.coef_[0])
                    coefs_alveolars.append(model_b.coef_[1])
                    coefs_velars.append(model_b.coef_[2])
            else:
                if multiway:
                    coefs_J.append(model_b.coef_J_)
                    coefs_K.append(model_b.coef_K_)
                else:
                    coefs.append(model_b.coef)
        if sparse:
            if problem == 'place':
                if multiway:
                    P_coefs_J_bilabials, kappa_J_bilabials, freq_J_bilabials = fleiss_kappa(coefs_J_bilabials, B)
                    P_coefs_J_alveolars, kappa_J_alveolars, freq_J_alveolars = fleiss_kappa(coefs_J_alveolars, B)
                    P_coefs_J_velars, kappa_J_velars, freq_J_velars = fleiss_kappa(coefs_J_velars, B)
                    P_coefs_K_bilabials, kappa_K_bilabials, freq_K_bilabials = fleiss_kappa(coefs_K_bilabials, B)
                    P_coefs_K_alveolars, kappa_K_alveolars, freq_K_alveolars = fleiss_kappa(coefs_K_alveolars, B)
                    P_coefs_K_velars, kappa_K_velars, freq_K_velars = fleiss_kappa(coefs_K_velars, B)

                    if write_to_file:
                        file_name = model.__class__.__name__ + '_' + problem + '_sparse_weights.csv'
                        with open(file_name, 'a') as file:
                            file.write(','.join([str(kappa_J_bilabials), str(kappa_K_bilabials),
                                                 *[str(p) for p in P_coefs_J_bilabials],
                                                 *[str(p) for p in P_coefs_K_bilabials],
                                                 *[str(f) for f in freq_J_bilabials],
                                                 *[str(f) for f in freq_K_bilabials],
                                                 ]))
                            file.write('\n')
                            file.write(','.join([str(kappa_J_alveolars), str(kappa_K_alveolars),
                                                 *[str(p) for p in P_coefs_J_alveolars],
                                                 *[str(p) for p in P_coefs_K_alveolars],
                                                 *[str(f) for f in freq_J_alveolars],
                                                 *[str(f) for f in freq_K_alveolars],
                                                 ]))
                            file.write('\n')
                            file.write(','.join([str(kappa_J_velars), str(kappa_K_velars),
                                                 *[str(p) for p in P_coefs_J_velars],
                                                 *[str(p) for p in P_coefs_K_velars],
                                                 *[str(f) for f in freq_J_velars],
                                                 *[str(f) for f in freq_K_velars],
                                                 ]))
                            file.write('\n')

                    return kappa_J_bilabials, P_coefs_J_bilabials, freq_J_bilabials, \
                           kappa_K_bilabials, P_coefs_K_bilabials, freq_K_bilabials, \
                           kappa_J_alveolars, P_coefs_J_alveolars, freq_J_alveolars, \
                           kappa_K_alveolars, P_coefs_K_alveolars, freq_K_alveolars, \
                           kappa_J_velars, P_coefs_J_velars, freq_J_velars, \
                           kappa_K_velars, P_coefs_K_velars, freq_K_velars

                P_coefs_bilabials, kappa_bilabials, freq_bilabials = fleiss_kappa(coefs_bilabials, B)
                P_coefs_alveolars, kappa_alveolars, freq_alveolars = fleiss_kappa(coefs_alveolars, B)
                P_coefs_velars, kappa_velars, freq_velars = fleiss_kappa(coefs_velars, B)

                if write_to_file:
                    file_name = model.__class__.__name__ + '_' + problem + '_sparse_weights.csv'
                    with open(file_name, 'a') as file:
                        file.write(','.join([str(kappa_bilabials), *[str(p) for p in P_coefs_bilabials],
                                             *[str(f) for f in freq_bilabials]]))
                        file.write('\n')
                        file.write(','.join([str(kappa_alveolars), *[str(p) for p in P_coefs_alveolars],
                                             *[str(f) for f in freq_alveolars]]))
                        file.write('\n')
                        file.write(','.join([str(kappa_velars), *[str(p) for p in P_coefs_velars],
                                             *[str(f) for f in freq_velars]]))
                        file.write('\n')

                return kappa_bilabials, P_coefs_bilabials, freq_bilabials, \
                            kappa_alveolars, P_coefs_alveolars, freq_alveolars, \
                            kappa_velars, P_coefs_velars, freq_velars

            if multiway:
                P_coefs_J, kappa_J, freq_J = fleiss_kappa(coefs_J, B)
                P_coefs_K, kappa_K, freq_K = fleiss_kappa(coefs_K, B)

                if write_to_file:
                    file_name = model.__class__.__name__ + '_' + problem + '_sparse_weights.csv'
                    with open(file_name, 'a') as file:
                        file.write(','.join(
                            [str(kappa_J), str(kappa_K), *[str(p) for p in P_coefs_J], *[str(p) for p in P_coefs_K],
                             *[str(f) for f in freq_J], *[str(f) for f in freq_K]]))
                        file.write('\n')
                return kappa_J, P_coefs_J, freq_J, kappa_K, P_coefs_K, freq_K

            P_coefs, kappa, freq = fleiss_kappa(coefs, B)
            if write_to_file:
                file_name = model.__class__.__name__ + '_' + problem + '_sparse_weights.csv'
                with open(file_name, 'a') as file:
                    file.write(','.join([str(kappa), *[str(p) for p in P_coefs], *[str(f) for f in freq]]))
                    file.write('\n')
            return kappa, P_coefs, freq

        else:
            if problem == 'place':
                if multiway:
                    mean_coefs_J_bilabials = np.mean(coefs_J_bilabials, axis=0)
                    var_coefs_J_bilabials = np.var(coefs_J_bilabials, axis=0)
                    mean_coefs_J_alveolars = np.mean(coefs_J_alveolars, axis=0)
                    var_coefs_J_alveolars = np.var(coefs_J_alveolars, axis=0)
                    mean_coefs_J_velars = np.mean(coefs_J_velars, axis=0)
                    var_coefs_J_velars = np.var(coefs_J_velars, axis=0)

                    mean_coefs_K_bilabials = np.mean(coefs_K_bilabials, axis=0)
                    var_coefs_K_bilabials = np.var(coefs_K_bilabials, axis=0)
                    mean_coefs_K_alveolars = np.mean(coefs_K_alveolars, axis=0)
                    var_coefs_K_alveolars = np.var(coefs_K_alveolars, axis=0)
                    mean_coefs_K_velars = np.mean(coefs_K_velars, axis=0)
                    var_coefs_K_velars = np.var(coefs_K_velars, axis=0)

                    if write_to_file:
                        file_name = model.__class__.__name__ + '_' + problem + '_non_sparse_weights.csv'
                        with open(file_name, 'a') as file:
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_J_bilabials], *[str(c) for c in var_coefs_J_bilabials]]))
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_K_bilabials], *[str(c) for c in var_coefs_K_bilabials]]))
                            file.write('\n')
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_J_alveolars], *[str(c) for c in var_coefs_J_alveolars]]))
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_K_alveolars], *[str(c) for c in var_coefs_K_alveolars]]))
                            file.write('\n')
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_J_velars], *[str(c) for c in var_coefs_J_velars]]))
                            file.write(','.join(
                                [*[str(c) for c in mean_coefs_K_velars], *[str(c) for c in var_coefs_K_velars]]))
                            file.write('\n')
                    return mean_coefs_J_bilabials, var_coefs_J_bilabials, mean_coefs_K_bilabials, var_coefs_K_bilabials, \
                           mean_coefs_J_alveolars, var_coefs_J_alveolars, mean_coefs_K_alveolars, var_coefs_K_alveolars, \
                           mean_coefs_J_velars, var_coefs_J_velars, mean_coefs_K_velars, var_coefs_K_velars

                mean_coefs_bilabials = np.mean(coefs_bilabials, axis=0)
                var_coefs_bilabials = np.var(coefs_bilabials, axis=0)
                mean_coefs_alveolars = np.mean(coefs_alveolars, axis=0)
                var_coefs_alveolars = np.var(coefs_alveolars, axis=0)
                mean_coefs_velars = np.mean(coefs_velars, axis=0)
                var_coefs_velars = np.var(coefs_velars, axis=0)

                if write_to_file:
                    file_name = model.__class__.__name__ + '_' + problem + '_non_sparse_weights.csv'
                    with open(file_name, 'a') as file:
                        file.write(
                            ','.join([*[str(c) for c in mean_coefs_bilabials], *[str(c) for c in var_coefs_bilabials]]))
                        file.write('\n')
                        file.write(
                            ','.join([*[str(c) for c in mean_coefs_alveolars], *[str(c) for c in var_coefs_alveolars]]))
                        file.write('\n')
                        file.write(
                            ','.join([*[str(c) for c in mean_coefs_velars], *[str(c) for c in var_coefs_velars]]))
                        file.write('\n')
                return mean_coefs_bilabials, var_coefs_bilabials, mean_coefs_alveolars, var_coefs_alveolars, mean_coefs_velars, var_coefs_velars

            if multiway:
                mean_coefs_J = np.mean(coefs_J, axis=0)
                var_coefs_J = np.var(coefs_J, axis=0)
                mean_coefs_K = np.mean(coefs_K, axis=0)
                var_coefs_K = np.var(coefs_K, axis=0)

                if write_to_file:
                    file_name = model.__class__.__name__ + '_' + problem + '_non_sparse_weights.csv'
                    with open(file_name, 'a') as file:
                        file.write(','.join([*[str(c) for c in mean_coefs_J], *[str(c) for c in var_coefs_J]]))
                        file.write(','.join([*[str(c) for c in mean_coefs_K], *[str(c) for c in var_coefs_K]]))
                        file.write('\n')
                return mean_coefs_J, var_coefs_J, mean_coefs_K, var_coefs_K

            mean_coefs = np.mean(coefs, axis=0)
            var_coefs = np.var(coefs, axis=0)
            if write_to_file:
                file_name = model.__class__.__name__ + '_' + problem + '_non_sparse_weights.csv'
                with open(file_name, 'a') as file:
                    file.write(','.join([*[str(c) for c in mean_coefs], *[str(c) for c in var_coefs]]))
                    file.write('\n')

            return mean_coefs, var_coefs
