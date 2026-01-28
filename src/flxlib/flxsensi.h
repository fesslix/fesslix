/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2026 Wolfgang Betz
 *
 * Fesslix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fesslix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fesslix.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "flxVec.h"
#include "flxrandom_base.h"


class flx_sensi_batch {
    private:
        std::valarray<vdouble> x_record;
        vdouble y_record;
    public:
        flx_sensi_batch(const tuint x_dim);
        void record_value(const flxVec& x_vec, const tdouble y);

        vdouble& get_y_ref() { return y_record; }

};

class flx_sensi_splitter_el {
    private:
        size_t Nbatches;        // there is always one batch reserved for NaN (the last one)
        std::vector<tdouble> splits;

    public:
        flx_sensi_splitter_el(const size_t Nsplit, const tdouble* vvecp, const size_t n);
        const size_t get_batch_index(const tdouble val) const;
        const size_t get_Nbatches() const { return Nbatches; }
};

class flx_sensi_splitter {
    private:
        size_t Nbatches;
        const tuint x_dim;
        std::valarray<flx_sensi_splitter_el*> splits;
        std::valarray<flx_sensi_batch>* batches;
    public:
        flx_sensi_splitter(const size_t Nsplit, const tuint x_dim, const std::valarray<flxVec*>& vvecp, const size_t n);
        ~flx_sensi_splitter();
        void record_value(const flxVec& x_vec, const tdouble y);
        const tdouble eval(const pdouble &y_mean, const pdouble &y_var);
        void eval_dist(flxVec& svec, FlxRndCreator &rndCreator, vdouble& y_record);
        const size_t get_Nbatches() const { return Nbatches; }
};


class flx_sensi_s1o {
    private:
        const std::string name;
        const size_t N_learn;
        const tuint x_dim;
        std::valarray<flxVec*> x_learn;
        flxVec* y_learn;
        tuint Nbrv;
        flx_sensi_splitter** brecord_vec;

        std::valarray<vdouble> x_record;
        vdouble y_record;

        tdouble eval_last;
        bool has_eval;

        void allocate_brecord();
    public:
        flx_sensi_s1o(const std::string name, const size_t N_learn, const tuint x_dim);
        virtual ~flx_sensi_s1o();

        void record_value(const flxVec& x_vec, const tdouble y);
        const tdouble eval();
        void eval_dist(flxVec& svec, FlxRndCreator &rndCreator);
        const tuint get_x_dim() const { return x_dim; }
};



