/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2025 Wolfgang Betz
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

#include <cmath>

#include "flxobjrandom.h"
#include "flxfunction_ope_calc.h"
#include "flxfunction_fun_calc.h"
#include "flxmtxfun_fun_calc.h"
#include "flxobjmtx.h"
#include "flxBayUp_obj.h"
#include "flxsensi.h"

#include <fstream>

flxBayUp* FlxObjReadSuS::lastSuS = nullptr;
FlxVoidBox<flx_sensi_s1o> sensi_s1o_box;

FlxRndCreator* flxPyRV::RndCreator_ptr = nullptr;

void FlxCreateObjReaders_RND::createObjReaders(FlxObjectReadBox* objReadBox) {
  objReadBox->insert("mci", new FlxObjReadMCI());
  objReadBox->insert("ips", new FlxObjReadIpS());
  objReadBox->insert("line_smpl", new FlxObjReadLineSmpl());
  objReadBox->insert("sus", new FlxObjReadSuS());
  objReadBox->insert("form", new FlxObjReadFORM());
  objReadBox->insert("form_pdf", new FlxObjReadFORM_pdf());
  objReadBox->insert("partial_derivative", new FlxObjReadFORM(true));
  objReadBox->insert("kde", new FlxObjReadKDE());
  objReadBox->insert("mcs_sensitivities", new FlxObjReadMCSsensi());
  objReadBox->insert("statsmp", new FlxObjReadStatSmp());
  objReadBox->insert("sortsmp", new FlxObjReadSortSmp());
  objReadBox->insert("smpplot", new FlxObjReadSmpPlot());
  objReadBox->insert("qq_plot", new FlxObjReadQQplot());
  objReadBox->insert("form_beta_sensitivities", new FlxObjReadFORMbetaSensitivities());
  objReadBox->insert("sus_level_info", new FlxObjReadSus_level_info());
  objReadBox->insert("sensi_s1o_new", new FlxObjReadSensi_s1o_new());
  objReadBox->insert("sensi_s1o_add", new FlxObjReadSensi_s1o_add());
  objReadBox->insert("sensi_s1o_dist", new FlxObjReadSensi_s1o_dist());
}

void FlxCreateObjReaders_RND::createFunReaders(FlxData* dataBox)
{
  dataBox->ConstantBox.declareC("sss_iter");
  
  dataBox->FunBox.insert("cdf_smp", new FunReadFunSmpCDF() );
  dataBox->FunBox.insert("sensi_s1o_eval", new FunReadFunSensi_s1o_eval() );
}


// #################################################################################
// random variables
// #################################################################################

RBRV_entry_RV_base * parse_py_obj_as_rv(py::dict config, const bool name_required, const tuint iID, const std::string family, std::string descr)
{
    RBRV_entry_RV_base* rv_ptr = nullptr;
    // retrieve type
        const std::string rv_type = parse_str_as_word(parse_py_para_as_string("type",config,true),true);
    // retrieve name
        std::string rv_name = family + parse_str_as_word(parse_py_para_as_string("name",config,name_required,"name_unspecified"),true);
    // select rbrv-class based on type
        if (rv_type=="stdn") {
            rv_ptr = new RBRV_entry_RV_stdN(rv_name,iID,config);
        } else if (rv_type=="normal") {
            rv_ptr = new RBRV_entry_RV_normal(rv_name,iID,config);
        } else if (rv_type=="logn") {
            rv_ptr = new RBRV_entry_RV_lognormal(rv_name,iID,config);
        } else if (rv_type=="uniform") {
            rv_ptr = new RBRV_entry_RV_uniform(rv_name,iID,config);
        // } else if (rv_type=="gumbel") {
        //     rv_ptr = new RBRV_entry_read_Gumbel(rv_name,iID,config);
        // } else if (rv_type=="normal_trunc") {
        //     rv_ptr = new RBRV_entry_read_normal_trunc(rv_name,iID,config);
        } else if (rv_type=="beta") {
            rv_ptr = new RBRV_entry_RV_beta(rv_name,iID,config);
        // } else if (rv_type=="exponential") {
        //     rv_ptr = new RBRV_entry_read_exponential(rv_name,iID,config);
        // } else if (rv_type=="gamma") {
        //     rv_ptr = new RBRV_entry_read_gamma(rv_name,iID,config);
        // } else if (rv_type=="poisson") {
        //     rv_ptr = new RBRV_entry_read_Poisson(rv_name,iID,config);
        // } else if (rv_type=="binomial") {
        //     rv_ptr = new RBRV_entry_read_Binomial(rv_name,iID,config);
        // } else if (rv_type=="cauchy") {
        //     rv_ptr = new RBRV_entry_read_Cauchy(rv_name,iID,config);
        // } else if (rv_type=="weibull") {
        //     rv_ptr = new RBRV_entry_read_Weibull(rv_name,iID,config);
        // } else if (rv_type=="chisquared") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(true,rv_name,iID,config);
        // } else if (rv_type=="chi") {
        //     rv_ptr = new RBRV_entry_read_ChiSquared(false,rv_name,iID,config);
        } else if (rv_type=="studentst") {
            rv_ptr = new RBRV_entry_RV_StudentsT(rv_name,iID,config);
        } else if (rv_type=="studentstgen") {
            rv_ptr = new RBRV_entry_RV_StudentsT_generalized(rv_name,iID,config);
        } else if (rv_type=="logt") {
            rv_ptr = new RBRV_entry_RV_logt(rv_name,iID,config);
        // } else if (rv_type=="laplace") {
        //     rv_ptr = new RBRV_entry_read_Laplace(rv_name,iID,config);
        // } else if (rv_type=="usertransform") {
        //     rv_ptr = new RBRV_entry_read_UserTransform(rv_name,iID,config);
        // } else if (rv_type=="truncated") {
        //     rv_ptr = new RBRV_entry_read_Truncated(rv_name,iID,config);
        // } else if (rv_type=="maxmintransform") {
        //     rv_ptr = new RBRV_entry_read_maxminTransform(rv_name,iID,config);
        } else if (rv_type=="quantiles") {
            rv_ptr = new RBRV_entry_RV_quantiles(rv_name,iID,config);
        // } else if (rv_type=="bayda") {
        //     rv_ptr = new RBRV_entry_read_bayDA(rv_name,iID,config);
        } else {
            std::ostringstream ssV;
            ssV << "Unknown random variable type '" << rv_type << "'.";
            throw FlxException("flxPyRV::flxPyRV_50", ssV.str() );
        }
    return rv_ptr;
}


flxPyRV::flxPyRV(py::dict config)
: rv_ptr(nullptr), mem_managed(true)
{
    rv_ptr_ = parse_py_obj_as_rv(config, false, 0, "", "flx.rv");
    rv_ptr = rv_ptr_;
}

flxPyRV::flxPyRV(RBRV_entry* rv_ptr, const bool mem_managed)
: rv_ptr(rv_ptr), rv_ptr_(dynamic_cast<RBRV_entry_RV_base*>(rv_ptr)), mem_managed(mem_managed)
{

}

flxPyRV::flxPyRV(flxPyRV& rhs)
: rv_ptr(rhs.rv_ptr), rv_ptr_(rhs.rv_ptr_), mem_managed(rhs.mem_managed)
{
    if (mem_managed) {
        throw FlxException_NotImplemented("flxPyRV::flxPyRV(flxPyRV& rhs)");
    }
}

flxPyRV::flxPyRV(flxPyRV&& rhs)
: rv_ptr(rhs.rv_ptr), rv_ptr_(rhs.rv_ptr_), mem_managed(rhs.mem_managed)
{
    rhs.rv_ptr = nullptr;
    rhs.rv_ptr_ = nullptr;
    rhs.mem_managed = false;
}

flxPyRV::~flxPyRV()
{
    if (rv_ptr) {
        if (mem_managed) {
            delete rv_ptr;
        }
    }
}

void flxPyRV::ensure_is_a_basic_rv()
{
    if (rv_ptr_==nullptr) {
        throw FlxException_NeglectInInteractive("flxPyRV::ensure_is_a_basic_rv", "'" + rv_ptr->name + "' is not a basic random variable.");
    }
}

const std::string flxPyRV::get_name() const
{
    return rv_ptr->name;
}

const std::string flxPyRV::get_type() const
{
    return rv_ptr->get_type();
}

const tdouble flxPyRV::get_value() const
{
    return rv_ptr->get_value();
}

const tdouble flxPyRV::x2y(const tdouble x_val)
{
    rv_ptr->eval_para();
    return rv_ptr->transform_x2y(x_val);
}

const tdouble flxPyRV::y2x(const tdouble y_val)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    return rv_ptr_->transform_y2x(y_val);
}

const tdouble flxPyRV::sample()
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    const tdouble y = get_RndCreator().gen_smp();
    return rv_ptr_->transform_y2x(y);
}

void flxPyRV::sample_array(py::array_t<tdouble> arr)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();

    // Access the data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* res_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Generate the samples (in standard Normal space)
    flxVec res_vec(res_ptr,size,false,false);
    get_RndCreator().gen_smp(res_vec);

    // transform the samples to original space
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr_->transform_y2x(res_ptr[i]);    // TODO avoid re-evaluating the parameters of the random variable
    }
}

const tdouble flxPyRV::pdf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_pdf_x(x_val,safeCalc);
}

py::array_t<tdouble> flxPyRV::pdf_array(py::array_t<tdouble> arr, const bool safeCalc)
{
    rv_ptr->eval_para();
    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Allocate memory for the return array
    auto res_buf = py::array_t<tdouble>(size);

    // Get the buffer info to access the underlying return data
    py::buffer_info res_buf_info = res_buf.request();
    tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);

    // Fill the array with values
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr->calc_pdf_x(input_ptr[i],safeCalc);    // TODO avoid re-evaluating the parameters of the random variable
    }

    // Return the array
    return res_buf;
}

const tdouble flxPyRV::pdf_log(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_pdf_x_log(x_val,safeCalc);
}

const tdouble flxPyRV::cdf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_cdf_x(x_val,safeCalc);
}

py::array_t<tdouble> flxPyRV::cdf_array(py::array_t<tdouble> arr, const bool safeCalc)
{
    rv_ptr->eval_para();

    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    // Allocate memory for the return array
    auto res_buf = py::array_t<tdouble>(size);

    // Get the buffer info to access the underlying return data
    py::buffer_info res_buf_info = res_buf.request();
    tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);

    // Fill the array with values
    for (size_t i = 0; i < size; ++i) {
        res_ptr[i] = rv_ptr->calc_cdf_x(input_ptr[i],safeCalc);    // TODO avoid re-evaluating the parameters of the random variable
    }

    // Return the array
    return res_buf;


}

const tdouble flxPyRV::icdf(const tdouble p)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    return rv_ptr_->calc_icdf_x(p);
}

const tdouble flxPyRV::sf(const tdouble x_val, const bool safeCalc)
{
    rv_ptr->eval_para();
    return rv_ptr->calc_sf_x(x_val,safeCalc);
}

const tdouble flxPyRV::entropy()
{
    rv_ptr->eval_para();
    return rv_ptr->calc_entropy();
}

const tdouble flxPyRV::mean()
{
    rv_ptr->eval_para();
    return rv_ptr->get_mean_current_config();
}

const tdouble flxPyRV::sd()
{
    rv_ptr->eval_para();
    return rv_ptr->get_sd_current_config();
}

const tdouble flxPyRV::median()
{
    rv_ptr->eval_para();
    return rv_ptr->get_median_current_config();
}

const tdouble flxPyRV::mode()
{
    rv_ptr->eval_para();
    return rv_ptr->get_mode_current_config();
}

const bool flxPyRV::check_x(const tdouble xV)
{
    rv_ptr->eval_para();
    return rv_ptr->check_x(xV);
}

const tdouble flxPyRV::get_HPD(const tdouble p)
{
    ensure_is_a_basic_rv();
    rv_ptr_->eval_para();
    return rv_ptr_->get_HPD(p);
}

py::dict flxPyRV::info()
{
    rv_ptr->eval_para();
    return rv_ptr->info();
}

FlxRndCreator& flxPyRV::get_RndCreator()
{
  if (RndCreator_ptr) {
    return *RndCreator_ptr;
  } else {
    throw FlxException("flxPyRV::get_RndCreator","Please start the engine.");
  }
}



// #################################################################################
// post-processors
// #################################################################################

post_proc_mean_double::post_proc_mean_double(const tuint colID)
: sum(ZERO), N(0), colID(colID)
{

}

void post_proc_mean_double::append_data(const flxVec& vec_full)
{
  sum += vec_full[colID];
  ++N;
}

py::dict post_proc_mean_double::eval()
{
    py::dict res;
    res["N"] = N;
    res["mean"] = sum/N;

    return res;
}

post_proc_mean_pdouble::post_proc_mean_pdouble(const tuint colID)
: N(0), colID(colID)
{

}

void post_proc_mean_pdouble::append_data(const flxVec& vec_full)
{
  sum += vec_full[colID];
  ++N;
}

py::dict post_proc_mean_pdouble::eval()
{
    py::dict res;
    res["N"] = N;
    res["mean"] = sum.cast2double()/N;

    return res;
}

post_proc_mean_qdouble::post_proc_mean_qdouble(const tuint colID, const size_t NpV,const bool ppb)
: sum(NpV,ppb), colID(colID)
{

}

void post_proc_mean_qdouble::append_data(const flxVec& vec_full)
{
  sum += vec_full[colID];
}

py::dict post_proc_mean_qdouble::eval()
{
    py::dict res;
    res["mean"] = sum.get_average();
    res["N"] = sum.get_N();
    return res;
}

post_proc_mean_vdouble::post_proc_mean_vdouble(const tuint colID)
: colID(colID)
{

}

void post_proc_mean_vdouble::append_data(const flxVec& vec_full)
{
  sum += vec_full[colID];
}

py::dict post_proc_mean_vdouble::eval()
{
    py::dict res;
    const tdouble mu = sum.get_mean();
    res["mean"] = mu;
    const tdouble var = sum.get_variance();
    res["var"] = var;
    const tdouble sd  = sqrt(var);
    res["sd"] = sd;
    const size_t N = sum.get_size();
    res["N"] = N;
    // assembe distribution of mean
      py::dict rv_mean_config;
      rv_mean_config["dof"] = tdouble(N)-ONE;
      rv_mean_config["loc"] = mu;
      rv_mean_config["scale"] = sd;
      rv_mean_config["eval_once"] = true;
      RBRV_entry_RV_StudentsT_generalized* rv_mean = new RBRV_entry_RV_StudentsT_generalized("unspecified",0,rv_mean_config);
      res["rv_mean"] = flxPyRV(rv_mean,true);
    return res;
}

post_proc_mean_reliability::post_proc_mean_reliability(const tuint colID)
: N(0), H(0), colID(colID)
{

}

void post_proc_mean_reliability::append_data(const flxVec& vec_full)
{
  ++N;
  if (vec_full[colID]<=ZERO) {
    ++H;
  }
}

py::dict post_proc_mean_reliability::eval()
{
    py::dict res;
    res["N"] = N;
    res["H"] = H;
    res["mean_freq"] = tdouble(H)/tdouble(N);
    // assembe distribution of pf
      py::dict rv_pf_config;
      rv_pf_config["alpha"] = tdouble(H)+ONE;
      rv_pf_config["beta"] = tdouble(N-H)+ONE;
      rv_pf_config["eval_once"] = true;
      RBRV_entry_RV_beta* rv_pf = new RBRV_entry_RV_beta("unspecified",0,rv_pf_config);
      res["rv_pf"] = flxPyRV(rv_pf,true);
    res["mean_bayes"] = rv_pf->get_mean_current_config();
    return res;
}

// #################################################################################
// dataBox
// #################################################################################

flxDataBox::flxDataBox(const tuint M_in, const tuint M_out)
: M_in(M_in), M_out(M_out), M(M_in+M_out),
    vec_full(M_in+M_out),
    vec_out(vec_full.get_tmp_vptr(),M_out),
    vec_in(vec_full.get_tmp_vptr()+M_out,M_in),
    fstream(nullptr), fs_N_col(0), fs_cols(nullptr), fs_binary(true),
    mem_N_reserved(0), mem_N(0), mem_ptr(nullptr), mem_cols(nullptr)
{

}

flxDataBox::~flxDataBox()
{
  free_mem();
  close_file();
  for (post_proc_base* pp_ptr : pp_vec) {
      if (pp_ptr) delete pp_ptr;
  }
}

void flxDataBox::ensure_M(const tuint M_) const
{
  if (M_!=M) {
    throw FlxException("flxDataBox::ensure_M", "Mismatch in dimension of input vector.");
  }
}

void flxDataBox::ensure_M_in(const tuint M_) const
{
  if (M_!=M_in) {
    throw FlxException("flxDataBox::ensure_M_in", "Mismatch in dimension of input vector.");
  }
}

void flxDataBox::ensure_M_out(const tuint M_) const
{
  if (M_!=M_out) {
    throw FlxException("flxDataBox::ensure_M_out", "Mismatch in dimension of output vector.");
  }
}

void flxDataBox::append_data()
{
  // write to file
    if (fstream) {
      if (fs_binary) {
        for (tuint i=0;i<fs_N_col;++i) {
          const tfloat value = vec_full[fs_cols[i]];
          fstream->write(reinterpret_cast<const char*>(&value), sizeof(tfloat));
        }
      } else {
        for (tuint i=0;i<fs_N_col;++i) {
          if (i>0) (*fstream) << ", ";
          (*fstream) << GlobalVar.Double2String(vec_full[fs_cols[i]]);
        }
        (*fstream) << std::endl;
      }
    }

  // write to memory
    if (mem_ptr) {
      if (mem_N>=mem_N_reserved) {
        throw FlxException("flxDataBox::append_data_20", "Memory of dataBox is full.");
      }
      tfloat* mem_ptr_pos = mem_ptr + mem_N;
      for (tuint i=0;i<mem_N_col;++i) {
        (*mem_ptr_pos) = vec_full[mem_cols[i]];
        mem_ptr_pos += mem_N_reserved;
      }
      ++mem_N;
    }

  // handle post-processors
    for (post_proc_base* pp_ptr : pp_vec) {
        pp_ptr->append_data(vec_full);
    }
}

tuint * flxDataBox::process_col_input(tuint& N_col, py::dict config)
{
  tuint* col_ptr = nullptr;
  try {
    // assign columns to write
      std::string col_str = "all";
      if (config.contains("cols")) {
        py::object obj = config["cols"];
        if (py::isinstance<py::list>(obj)) {
          // columns provided as list
            py::list col_lst = obj.cast<py::list>();
            N_col = col_lst.size();
            col_ptr = new tuint[N_col];
            for (tuint i=0;i<N_col;++i) {
              const tuint col_id = col_lst[i].cast<tuint>();
              if (col_id>=M) {
                  std::ostringstream ssV;
                  ssV << "ID of column (" << col_id << ") must be smaller than " << M;
                  throw FlxException("flxDataBox::process_col_input_01", ssV.str());
              }
              col_ptr[i] = col_id;
            }
        } else {
          // columns provided as string
            col_str = parse_py_obj_as_string(obj,"flxDataBox::process_col_input_02");
        }
      }
      if (col_ptr==nullptr) {
        // process string
          if (col_str=="all") {
            N_col = M;
            col_ptr = new tuint[N_col];
            for (tuint i=0;i<N_col;++i) {
              col_ptr[i] = i;
            }
          } else if (col_str=="all_in") {
            N_col = M_in;
            col_ptr = new tuint[N_col];
            for (tuint i=0;i<N_col;++i) {
              col_ptr[i] = i+M_out;
            }
          } else if (col_str=="all_out") {
            N_col = M_out;
            col_ptr = new tuint[N_col];
            for (tuint i=0;i<N_col;++i) {
              col_ptr[i] = i;
            }
          } else {
             throw FlxException("flxDataBox::process_col_input_03", "Unknown keyword for 'col': " + col_str);
          }
      }
  } catch (FlxException& e) {
    if (col_ptr) {
      delete [] col_ptr;
      col_ptr = nullptr;
    }
    throw;
  }
  return col_ptr;
}

const tuint flxDataBox::extract_colID(py::object col)
{
  tuint colID = 0;
  // type: dict
    if (py::isinstance<py::dict>(col)) {
      py::dict dict = py::cast<py::dict>(col);
      const std::string set_str = parse_py_para_as_word("set",dict,true,true);
      const tuint id = parse_py_para_as_tuint("id",dict,true);
      if (set_str=="full") colID=id;
      else if (set_str=="out") colID=id;
      else if (set_str=="in") colID=id+M_out;
      else {
        throw FlxException("flxDataBox::extract_colID_01","Unknown keyword in 'set'");
      }
  // type: tuint
    } else {
      colID = parse_py_obj_as_tuint(col, "Value of key 'cols' in <dict> config.");
    }
  // check for consistency
    if (colID>=M) {
      throw FlxException("flxDataBox::extract_colID_99", "colID exceeds dimension of data-points.");
    }
    return colID;
}

const tuint flxDataBox::extract_colID_(py::dict config)
{
  if (config.contains("col")) {
    return extract_colID(config["col"]);
  } else {
    throw FlxException("flxDataBox::extract_colID_", "'col' not in <dict> 'config'.");
  }
}

void flxDataBox::write2mem(py::dict config)
{
  // checks
    if (mem_ptr) {
      throw FlxException("flxDataBox::write2mem_01","The dataBox is already linked to memory.");
    }
  // extract config from dict
    mem_N_reserved = parse_py_para_as_tulong("N_reserve",config,true);
  // get columns to store in memory
    mem_cols = process_col_input(mem_N_col, config);
  // allocate memory
    mem_ptr = new tfloat[mem_N_col*mem_N_reserved];
}

py::array_t<tfloat> flxDataBox::extract_col_from_mem(py::object col)
{
  const tuint colID = extract_colID(col);
  tfloat* mem_ptr_pos = mem_ptr + colID*mem_N_reserved;
  return py_wrap_array_no_ownership<tfloat>(mem_ptr_pos,mem_N);
}

void flxDataBox::free_mem()
{
  if (mem_ptr) {
    delete [] mem_ptr;
    mem_ptr = nullptr;
  }
  mem_N_reserved = 0;
  mem_N = 0;
  mem_N_col = 0;
  if (mem_cols) {
    delete [] mem_cols;
    mem_cols = nullptr;
  }
}

void flxDataBox::write2file(py::dict config)
{
  // checks
    if (fstream) {
      throw FlxException("flxDataBox::write2file_01","The dataBox is already linked to an output stream.");
    }
  // extract config from dict
    const std::string fname = parse_py_para_as_string("fname",config,true);
    const bool append = parse_py_para_as_bool("append",config,false,true);
    fs_binary = parse_py_para_as_bool("binary",config,false,true);
  // open file
    if ( append ) {
      fstream = new std::ofstream(fname.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::app );
    } else {
      fstream = new std::ofstream(fname.c_str(), std::ios_base::out | std::ios_base::binary );
    }
  try {
    // ensure file is open
      if ( fstream == nullptr || ! fstream->is_open() ) {
        std::ostringstream ssV;
        ssV << "File (" << fname << ") could not be opened.";
        throw FlxException("flxDataBox::write2file_02", ssV.str() );
      }
    fs_cols = process_col_input(fs_N_col, config);
  } catch (FlxException& e) {
    // close file
    close_file();
    throw;
  }
}

void flxDataBox::read_from_file(py::dict config)
{
  // extract config from dict
    const std::string fname = parse_py_para_as_string("fname",config,true);
    const bool ifs_binary = parse_py_para_as_bool("binary",config,false,true);
  if (ifs_binary) {
    // open file
      std::ifstream ifs(fname.c_str(), std::ios_base::in | std::ios_base::binary );
    // ensure file is open
      if ( !ifs ) {
        std::ostringstream ssV;
        ssV << "File (" << fname << ") could not be opened.";
        throw FlxException("flxDataBox::read_from_file_01", ssV.str() );
      }
    // Determine file size
      ifs.seekg(0, std::ios::end);
      std::streamsize size = ifs.tellg();
      ifs.seekg(0, std::ios::beg);
    // Compute how many floats the file holds
      const size_t num_floats = size / sizeof(tfloat);
      if (num_floats % M != 0) {
        std::ostringstream ssV;
        ssV << "Total number of values in the file is not a multiple of M (" << M << ").";
        throw FlxException("flxDataBox::read_from_file_02a", ssV.str() );
      }
    // allocate memory for importing float values
      std::vector<tfloat> buffer(M);
    // import the entire file
      while (ifs.read(reinterpret_cast<char*>(buffer.data()), M * sizeof(float))) {
        // cast to tdouble
          for (tuint i=0;i<M;++i) {
            vec_full[i] = buffer[i];
          }
        // import data
          append_data();
      }

    // Check if there was an incomplete final chunk (shouldn't happen if total is multiple of M)
      if (ifs.gcount() > 0) {
        throw FlxException_Crude("flxDataBox::read_from_file_03a");
      }
  } else {
    ReadStream rs(fname.c_str());
    while (rs.check_eof()==false) {
      // read next data point
        for (tuint i = 0; i < M; ++i) {
          if (rs.check_eof()) {
            std::ostringstream ssV;
            ssV << "Total number of values in the file is not a multiple of M (" << M << ").";
            throw FlxException("flxDataBox::read_from_file_02b", ssV.str() );
          }
          vec_full[i] = rs.get_Double();
          if (!(rs.check_eof())) {
            // set next
              const char c = rs.whatIsNextChar();
              if (c == ',' || c==';') {
                rs.getChar();
              }
          }
        }
      // import data
        append_data();
    }
  }
}

void flxDataBox::close_file()
{
    if (fstream) {
      delete fstream;
      fstream = nullptr;
    }
    fs_N_col = 0;
    if (fs_cols) {
      delete [] fs_cols;
      fs_cols = nullptr;
    }
}

post_proc_base& flxDataBox::register_post_processor(py::dict config)
{
  // extract type from config
    const std::string pp_type = parse_py_para_as_word("type",config,true,true);
  // define the post-processor
    post_proc_base* pp_ptr = nullptr;
    try {
      if (pp_type=="mean_double") {
        const tuint colID = extract_colID_(config);
        pp_ptr = new post_proc_mean_double(colID);
      } else if (pp_type=="mean_pdouble") {
        const tuint colID = extract_colID_(config);
        pp_ptr = new post_proc_mean_pdouble(colID);
      } else if (pp_type=="mean_qdouble") {
        const tuint colID = extract_colID_(config);
        const tuint NpV = parse_py_para_as_tuintNo0("NpV",config,false,10000000);
        const bool ppb = parse_py_para_as_bool("bbp",config,false,false);
        pp_ptr = new post_proc_mean_qdouble(colID,NpV,ppb);
      } else if (pp_type=="vdouble") {
        const tuint colID = extract_colID_(config);
        pp_ptr = new post_proc_mean_vdouble(colID);
      } else if (pp_type=="reliability") {
        const tuint colID = extract_colID_(config);
        pp_ptr = new post_proc_mean_reliability(colID);
      } else {
        throw FlxException("flxDataBox::register_post_processor_99","Unknown type ('" + pp_type + "' for post-processor.");
      }
    } catch (FlxException& e) {
      if (pp_ptr) delete pp_ptr;
      throw;
    }
  // append it to pp_vec
    pp_vec.push_back(pp_ptr);
  return *pp_ptr;
}



//======================== Monte Carlo Integration ===========================

void flx_perform_MCS(const tulong N, FlxMtxFun_base& model, RBRV_constructor& sampler, flxDataBox& dbox)
{
  flxVec& vec_out = dbox.vec_out;
  flxVec& vec_in = dbox.vec_in;
  tdouble* ptr_in = vec_in.get_tmp_vptr();
  // consistency checks
    dbox.ensure_M_in(sampler.get_NOX());
    dbox.ensure_M_out(model.get_N());
  // perform the Monte Carlo simulation
    for (tulong i=0;i<N;++i) {
      sampler.gen_smp();
      model.eval();
      vec_out = model.get_res_vec();
      sampler.get_x_Vec(ptr_in);
      dbox.append_data();
    }
}


void FlxObjMCI::FirstThingsFirst(RBRV_constructor& RndBox)
{
  GlobalVar.slogcout(4) << "mci: performing a Monte Carlo integration. (N=" << GlobalVar.Double2String(Np) << ")" << std::endl; 
}

void FlxObjMCI::log_AddResInfo(std::ostream& sout, const tdouble hits, const tdouble Npd)
{
  output_Bayesian_credible_reliablity(sout, hits, Npd);
}

void FlxObjMCI::output_Bayesian_credible_reliablity(std::ostream& sout, const tdouble hits, const tdouble Npd)
{
  if (reliability) {
    // Bayesian estimation
      const tdouble alpha_prior = ONE;
      const tdouble beta_prior = ONE;
      RBRV_entry_RV_beta betad("p",0,false,new FlxFunction(new FunNumber(hits+alpha_prior)),new FlxFunction(new FunNumber(Npd-hits+beta_prior)),NULL,NULL,true); 
      betad.eval_para();
      sout << "  mean(Pf) (Bayesian):      " << GlobalVar.Double2String(betad.get_mean_current_config()) << std::endl;
      sout << "  C.o.V.(Pf) (Bayesian): " << GlobalVar.Double2String(betad.get_sd_current_config()/betad.get_mean_current_config()) << std::endl;
      sout << "  Estim. C.o.V. (MLE):   " << GlobalVar.Double2String(sqrt((ONE-theResult)/(theResult*Npd))) << std::endl;
      // update _bne-const
        std::string res_name = data->ConstantBox.get(&theResult);
        res_name += "_bne";
        tdouble* theResult_bne = data->ConstantBox.get(res_name,true);
        *theResult_bne = betad.get_mean_current_config();
    // output credible intervals
    FlxSMtx* pcM = data->ConstMtxBox.get(pc->eval(),true);
    if (pcM->get_nrows()*pcM->get_ncols()>0) {
      // where to store the credible intervals (two-sided) to?:
        tuint Nentries = pcM->get_nrows()*pcM->get_ncols();
        tuint Ncols = 3;
        tdouble *crp = data->ConstMtxBox.get_Mtx("mci_credible_2sided",Nentries,Ncols);
        tuint crp_id = 0;
      sout << "  Credible intervals (two-sided) of probability of failure:" << std::endl;
      sout << "                          Bayesian                Normal approx." << std::endl;
      for ( tuint i=0; i<pcM->get_nrows(); ++i ) {
        for ( tuint j=0; j<pcM->get_ncols(); ++j ) {
          const tdouble pcT = pcM->operator()(i,j);
          // Calculate the credible interval
            sout << "       " << GlobalVar.Double2String(pcT,false,4,6) << ":    ";
            crp[crp_id++] = pcT;
            tdouble kc = ONE-(ONE-pcT)/2;
            const tdouble bay_low = betad.Inv_cdf_x((ONE-pcT)/2);
            const tdouble bay_up = betad.Inv_cdf_x(kc);
            crp[crp_id++] = bay_low;
            crp[crp_id++] = bay_up;
            sout << "[" << format("%10.3e") % bay_low << ";" << format("%10.3e") % bay_up << "]";
            kc = rv_InvPhi(kc);
            tdouble lb = theResult - kc * sqrt(theResult*(ONE-theResult)/Npd);
            if (lb<ZERO) lb=ZERO;
            tdouble ub = theResult + kc * sqrt(theResult*(ONE-theResult)/Npd);
            if (ub>ONE) ub=ONE;
            sout << "    [" << format("%10.3e") % lb << ";" << format("%10.3e") % ub << "]" << std::endl;
        }
      }
      // where to store the credible intervals (upper bound) to?:
        Ncols = 2;
        crp = data->ConstMtxBox.get_Mtx("mci_credible_ubound",Nentries,Ncols);
        crp_id = 0;
      sout << "  Credible intervals (upper bound) of probability of failure:" << std::endl;
      sout << "                          Bayesian                Normal approx." << std::endl;
      for ( tuint i=0; i<pcM->get_nrows(); ++i ) {
        for ( tuint j=0; j<pcM->get_ncols(); ++j ) {
          const tdouble pcT = pcM->operator()(i,j);
          // Calculate the cerdible interval
            sout << "       " << GlobalVar.Double2String(pcT,false,4,6) << ":    ";
            crp[crp_id++] = pcT;
            const tdouble bay_ub = betad.Inv_cdf_x(pcT);
            crp[crp_id++] = bay_ub;
            sout << "       " << format("%10.3e") % bay_ub;
            const tdouble kc = rv_InvPhi(pcT);
            tdouble ub = theResult + kc * sqrt(theResult*(ONE-theResult)/Npd);
            if (ub>ONE) ub=ONE;
            sout << "                " << format("%10.3e") % ub << std::endl;
        }
      }
    } // end if: output credible intervals
  } // end if: reliability
}


void FlxObjMCI::task()
{
  // generate the random_creator
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrvsets->eval(true));
    RBRV_constructor RndBox(set_str_vec,data->rbrv_box);
  Np = funNp->cast2tulong();
  FirstThingsFirst(RndBox);
  #if FLX_KAHAN_MCI
    pdouble Integ, IntegI;
  #else
    tdouble Integ, IntegI;
  #endif
  tdouble hits, hitsI;
  // split the interval
    tulong NpI = Np;
    if (interv) {
      NpI = tulong(sqrt(tdouble(Np)));
    };
  // Perform the integration
    Integ=ZERO; hits=ZERO;
    // load the progress bar
      FlxProgress prg(*GlobalVar.get_cout(),!NOTdolog);
      prg.start(NpI);
    if (NpI == Np) {
      for (tulong i = 0; i < Np; ++i) {
        Integrationstep(Integ,hits,RndBox);
        prg.tick(i);
      }
    } else {
      for (tulong j = 0; j < NpI-1; ++j) {
        IntegI = ZERO; hitsI = ZERO;
        for (tulong i = 0; i < NpI; ++i) {
          Integrationstep(IntegI,hitsI,RndBox);
        }
        prg.tick(j);
        Integ+=IntegI;
        hits+=hitsI;
      }
      IntegI = ZERO; hitsI = ZERO;
      for (tulong i = 0; i < Np-(NpI-1)*NpI; ++i) {
        Integrationstep(IntegI,hitsI,RndBox);
      }
      Integ+=IntegI;
      hits+=hitsI;
    }
    // stop the progress-bar
      prg.stop();
    Integ/=tdouble(Np); 
    #if FLX_KAHAN_MCI
      theResult = Integ.cast2double();
    #else
      theResult = Integ;
    #endif
    GlobalVar.slogcout(4) << " Result of the Integration: " << GlobalVar.Double2String(theResult);
    if (reliability) {
      GlobalVar.slogcout(4) << " (" << hits << " hits)";
    }
    GlobalVar.slogcout(4) << std::endl;
    if (reliability && (theResult > ONE || theResult < ZERO) ) {
      GlobalVar.alert.alert("FlxObjMCI::task","Result does not seem to be a probability !!!");
    }
    log_AddResInfo(GlobalVar.slogcout(4),hits,tdouble(Np));
    LastThingsLast();
}

FlxObjMCI::~FlxObjMCI()
{
  delete funNp;
  delete fung;
  delete pc;
  delete rbrvsets;
}

FlxObjReadMCI::FlxObjReadMCI()
{
  AllDefParaBox->insert(new FlxOptionalParaFlxString("nataf","sim::rbrvsets",true));
  ParaBox.insert("rbrvsets", "sim::rbrvsets" );
  
  // optional parameters
    // interval split
      AllDefParaBox->insert(new FlxOptionalParaBool(true,"mci::interv"));
      ParaBox.insert("interv", "mci::interv" );
    // determining probability of failure?
      AllDefParaBox->insert(new FlxOptionalParaBool(false,"mci::reliability"));
      ParaBox.insert("reliability", "mci::reliability" );
    // credible intervals
      std::vector<FlxFunction*> vecV;
      vecV.push_back(new FlxFunction(new FunNumber(0.90)));
      vecV.push_back(new FlxFunction(new FunNumber(0.95)));
      FlxObjBase* blockV = new FlxObjMtxConstNewU(false, new FlxMtxConstFun("internal_mcirelpc"), vecV, 2, 1);
      FlxMtxConstFun* mtxvalue = new FlxMtxConstFun( "internal_mcirelpc", blockV );
      AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"mci::reliability::pc"));
      ParaBox.insert("pc", "mci::reliability::pc" );
}

void FlxObjReadMCI::read_MCIblock(tdoublePtr& theResult, FlxFunctionPtr& funNp, FlxFunctionPtr& fung, bool errSerious)
{
  reader->getChar('(',errSerious);
  // result-constant
    const std::string cname = reader->getWord(true,errSerious);
    data->ConstantBox.declareC(cname);
    theResult = data->ConstantBox.get(cname);
    reader->getChar(';',errSerious);
  // function for number of points
    funNp = new FlxFunction(funReader,errSerious);
    reader->getChar(';',errSerious);
  // function to integrate
    fung = new FlxFunction(funReader,errSerious);
  reader->getChar(')',errSerious);
}

FlxObjBase* FlxObjReadMCI::read()
{
  tdouble* d1=NULL; FlxFunction* funNp=NULL; FlxFunction* fung=NULL;
  try {
    read_MCIblock(d1,funNp,fung,false);
    read_optionalPara(false);
    return new FlxObjMCI(get_doLog(), *d1, funNp, fung, get_optPara_bool("interv"),get_verboseLog(),get_optPara_bool("reliability"),get_optPara_FlxMtxFun("pc"),get_optPara_FlxString("rbrvsets"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadMCI::read",1);
    if (funNp) delete funNp;
    if (fung) delete fung;
    throw;
  }
}

void FlxObjIpS::FirstThingsFirst(RBRV_constructor& RndBox)
{
  GlobalVar.slogcout(4) << "ips: performing an Importance sampling integration (N=" << GlobalVar.Double2String(Np) << ")" << std::endl; 
  sspace = sspaceG->generate_SS(RndBox);
  GlobalVar.slogcout(4) << "  Sampling space: ";
  sspace->print_info(GlobalVar.slogcout(4), verboseLog);
  GlobalVar.slogcout(4) << std::endl;
  sum_weights_hit.clear();
  sum_weights_nohit.clear();
  sum_weights_total.clear();
}

void FlxObjIpS::log_AddResInfo(std::ostream& sout, const tdouble hits, const tdouble Npd)
{
//   if (sspace->rely_confidence()>ZERO) {
//     const tdouble f = sspace->rely_confidence();
//     FlxObjMCI::log_AddResInfo(sout, hits, Npd*f);
//   }
  // get scaling factor for equivalent Monte Carlo samples (in case of truncated importance sampling domain)
    const tdouble f = (sspace->rely_confidence()>ZERO)?(sspace->rely_confidence()):ONE;
  if (reliability) {
    const tdouble is_result = (theResult>ONE)?ONE:theResult;
    sout << "  Experimental output (for research only):   " << std::endl;
    sout << "    sum of weights: " << GlobalVar.Double2String(sum_weights_hit.get_sum()) << " :: " << GlobalVar.Double2String(sum_weights_nohit.get_sum()) << " :: " << GlobalVar.Double2String(sum_weights_total.get_sum()) << std::endl;
    sout << "    mean of weights: " << GlobalVar.Double2String(sum_weights_hit.get_mean()*f) << " :: " << GlobalVar.Double2String(sum_weights_nohit.get_mean()*f) << " :: " << GlobalVar.Double2String(sum_weights_total.get_mean()*f) << std::endl;
    sout << "    variance of weights: " << GlobalVar.Double2String(sum_weights_hit.get_variance()*pow2(f)) << " :: " << GlobalVar.Double2String(sum_weights_nohit.get_variance()*pow2(f)) << " :: " << GlobalVar.Double2String(sum_weights_total.get_variance()*pow2(f)) << std::endl;
    const tdouble N1 = pow(sum_weights_hit.get_sum(),2)/sum_weights_hit.get_sum_of_squares();
    const tdouble N2 = pow(sum_weights_nohit.get_sum(),2)/sum_weights_nohit.get_sum_of_squares();
    const tdouble N3 = pow(sum_weights_total.get_sum(),2)/sum_weights_total.get_sum_of_squares();
    sout << "    effective sample size: " << 
        GlobalVar.Double2String(N1) << " (" << GlobalVar.Double2String(hits) << ") " << " :: " << 
        GlobalVar.Double2String(N2) << " (" << GlobalVar.Double2String(Npd-hits) << ") "  << " :: " << 
        " total: " << GlobalVar.Double2String(N1+N2) << " (" << GlobalVar.Double2String(Npd) << ") :: " <<
        GlobalVar.Double2String(N3) << " (" << GlobalVar.Double2String(Npd) << ") "  << " :: " <<  std::endl;
    sout << "    m1  " << GlobalVar.Double2String(is_result*N3*f) << "  ::  " << GlobalVar.Double2String(N3*f) << std::endl;
    sout << "    m2  " << GlobalVar.Double2String(N1) << "  ::  " << GlobalVar.Double2String(N1/is_result) << std::endl;
    const tdouble m2a_a = (hits<2)?(is_result*N3*f):N1;
    const tdouble m2a_b = (hits<2)?(N3*f):(N1/is_result);
    sout << "    m2a " << GlobalVar.Double2String(m2a_a) << "  ::  " << GlobalVar.Double2String(m2a_b) << std::endl;
    const tdouble m3a_rf = (ONE+1/sqrt(N3)) * ((sum_weights_total.get_mean()*f<ONE)?(sum_weights_total.get_mean()*f):ONE);
    const tdouble m3a_a = ((hits<2)?(is_result*N3*f):N1)*m3a_rf;
    const tdouble m3a_b = ((hits<2)?(N3*f):(N1/is_result))*m3a_rf;
    sout << "    m3a " << GlobalVar.Double2String(m3a_a) << "  ::  " << GlobalVar.Double2String(m3a_b) << std::endl;
    tdouble m3b_rf = (sum_weights_total.get_variance()*pow2(f)<GlobalVar.TOL())?ONE:(rv_Phi((sum_weights_total.get_mean()*f-0.95)/(f*sqrt(sum_weights_total.get_variance()/sum_weights_total.get_size()))));
    const tdouble m3b_a = ((hits<2)?(is_result*N3*f):N1)*m3b_rf;
    const tdouble m3b_b = ((hits<2)?(N3*f):(N1/is_result))*m3b_rf;
    sout << "    m3b " << GlobalVar.Double2String(m3b_a) << "  ::  " << GlobalVar.Double2String(m3b_b) << "  ::  " << GlobalVar.Double2String(m3b_rf) << std::endl;
//     const tdouble m3c_rf1 = log((sum_weights_total.get_variance()*pow2(f)<GlobalVar.TOL())?1e20:(rv_phi((ONE-sum_weights_total.get_mean()*f)/(f*sqrt(sum_weights_total.get_variance()/sum_weights_total.get_size())))/(f*sqrt(sum_weights_total.get_variance()/sum_weights_total.get_size()))));
//     const tdouble m3c_rf = rv_Phi((m3c_rf1-2.)/0.1);
//     const tdouble m3c_a = ((hits<2)?(is_result*N3*f):N1)*m3c_rf;
//     const tdouble m3c_b = ((hits<2)?(N3*f):(N1/is_result))*m3c_rf;
//     sout << "    m3c " << GlobalVar.Double2String(m3c_a) << "  ::  " << GlobalVar.Double2String(m3c_b) << "  ::  " << GlobalVar.Double2String(m3c_rf) << "  ::  " << GlobalVar.Double2String(m3c_rf1) << std::endl;
//     const tdouble m3d_var = ((sum_weights_total.get_sum_of_squares()*pow2(f)-sum_weights_total.get_size())/(sum_weights_total.get_size()-1))/sum_weights_total.get_size();
//     const tdouble m3d_rf = ONE-rv_Phi(log(m3d_var));
//     const tdouble m3d_a = ((hits<2)?(is_result*N3*f):N1)*m3d_rf;
//     const tdouble m3d_b = ((hits<2)?(N3*f):(N1/is_result))*m3d_rf;
//     sout << "    m3d " << GlobalVar.Double2String(m3d_a) << "  ::  " << GlobalVar.Double2String(m3d_b) << "  ::  " << GlobalVar.Double2String(m3d_rf) << "  ::  " << GlobalVar.Double2String(m3d_var) << std::endl;
    const tdouble m4_Nt = pow2(Npd)/sum_weights_total.get_sum_of_squares();
    const tdouble m4_rf = pow(m4_Nt/Npd,0.1);
    const tdouble m4_a = ((hits<2)?(is_result*N3*f):N1)*m4_rf;
    const tdouble m4_b = ((hits<2)?(N3*f):(N1/is_result))*m4_rf;
    sout << "    m4 " << GlobalVar.Double2String(m4_a) << "  ::  " << GlobalVar.Double2String(m4_b) << "  ::  " << GlobalVar.Double2String(m4_rf) << "  ::  " << GlobalVar.Double2String(m4_Nt) << std::endl;
    // compute and output the Bayesian reliability bounds
//       output_Bayesian_credible_reliablity(sout, is_result*N3*f, N3*f);    // m1
//       output_Bayesian_credible_reliablity(sout, N1, N1/is_result );    // m2  
//       output_Bayesian_credible_reliablity(sout, m2a_a, m2a_b );    // m2a 
//       output_Bayesian_credible_reliablity(sout, m3a_a, m3a_b );    // m3a 
      output_Bayesian_credible_reliablity(sout, m4_a, m4_b );    // m3b 
  }
}

#if FLX_KAHAN_MCI
void FlxObjIpS::Integrationstep(pdouble& Integ, tdouble& hits, RBRV_constructor& RndBox)
#else
void FlxObjIpS::Integrationstep(tdouble& Integ, tdouble& hits, RBRV_constructor& RndBox)
#endif
{  
  tdouble foverh=ZERO; 
  sspace->transform_space(foverh);
  tdouble t = fung->calc()*foverh;
  Integ+=t;
  sum_weights_total += foverh;
  if (t>0) {
    ++hits;
    sum_weights_hit += foverh;
  } else {
    sum_weights_nohit += foverh;
  }
}

FlxObjBase* FlxObjReadIpS::read()
{
  tdouble* d1=NULL; FlxFunction* funNp=NULL; FlxFunction* fung=NULL;
  FlxRndSamplingSpace_Generator_base* sspaceG = NULL;
  try {
    read_MCIblock(d1,funNp,fung,false);
    // read sampling space
      reader->getChar('(',false);
      sspaceG = FlxRndSamplingSpace_Generator_base::createSS(FlxRndSamplingSpace_Generator_base::get_sst(reader->getWord(true,false),false),false);
      reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjIpS(get_doLog(), *d1, funNp, fung, get_optPara_bool("interv"),get_verboseLog(),get_optPara_bool("reliability"),get_optPara_FlxMtxFun("pc"),sspaceG,get_optPara_FlxString("rbrvsets"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadIpS::read",1);
    if (funNp) delete funNp;
    if (fung) delete fung;
    if (sspaceG) delete sspaceG;
    throw;
  }
}


const tdouble FlxObjLineSmpl::LSF_call(const tdouble c, const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, tulong& N_LSF_calls)
{
  rv_prop = rv_base;
  rv_prop.add(betaVec,c);
  RndBoxp->set_smp(rv_prop);
  ++N_LSF_calls;
  const tdouble res = fung->calc();
  if (extended_ls) {
    hist_push(c,res);
  }
  return res;
}

const tdouble FlxObjLineSmpl::perform_line_search(const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, const tdouble tol, const tuint iter_max, tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV)
{
  if (use_bisec) {
    return perform_line_search_bisec(rv_base,rv_prop,betaVec,tol,iter_max,N_LSF_calls,fdright,found,startV,endV);
  } else {
    return perform_line_search_rgfsi(rv_base,rv_prop,betaVec,tol,iter_max,N_LSF_calls,fdright,found,startV,endV);
  }
}

const tdouble FlxObjLineSmpl::perform_line_search_rgfsi(const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, const tdouble tol, const tuint iter_max, tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV)
{
  found = false;
  // evaluate the initial value
    tdouble c_end = endV;
    tdouble g_2 = LSF_call(c_end,rv_base,rv_prop,betaVec,N_LSF_calls);
  // evaluate the 2nd starting value
    tdouble c_start = startV;
    tdouble g_1 = LSF_call(c_start,rv_base,rv_prop,betaVec,N_LSF_calls);
  // start the iteration
    tuint c = 0;
    tdouble c_res=ZERO, g_3=ZERO;
    while (c<iter_max) {
      if (g_1*g_2<=ZERO) {        // Pegasus-algorithm
        // evaluate a new point
          c_res = (c_start*g_2-c_end*g_1)/(g_2-g_1);
          g_3 = LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls);
        if ( g_2*g_3<ZERO ) {
          c_start = c_end;
          g_1 = g_2;
          c_end = c_res;
          g_2 = g_3;
        } else {
          const tdouble m = g_2/(g_2+g_3);
          g_1 = m*g_1;
          c_end = c_res;
          g_2 = g_3;
        }
      } else {                        // Sekantenverfahren
        c_res = c_end - ((c_end-c_start)/(g_2-g_1))*g_2;
        if (fabs(c_res)>tdouble(50)) {
          found = false;
          return c_res;
        }
        g_3 = LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls);
        g_1 = g_2;
        c_start = c_end;
        g_2 = g_3;
        c_end = c_res;
      }
      // check for convergence
        if (fabs(g_3)<=tol) {
          found = true;
          break;
        }
        if (fabs(c_end-c_start)<=tol) break;
      ++c;
    }
  // try to deduce information about the topology
    fdright = ( g_1>g_2 );
    if (c_start>c_end) { fdright = !fdright; }
  // maximum number of line-search iterations reached?
    if (c>=iter_max) {
      if (verboseLog) {
        GlobalVar.alert.alert("FlxObjLineSmpl::perform_line_search_rgfsi","Maximum number of line-search iterations reached.");
      }
    }
  if (fabs(g_3)<=tol) return c_res;
  if (g_1*g_2<=ZERO && c<iter_max) found=true;
  return (c_start+c_end)/2;
}

const tdouble FlxObjLineSmpl::perform_line_search_bisec(const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, const tdouble tol, const tuint iter_max, tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV)
{
  found = false;
  // evaluate the initial value
    tdouble c_end = endV;
    tdouble g_2 = LSF_call(c_end,rv_base,rv_prop,betaVec,N_LSF_calls);
  // evaluate the 2nd starting value
    tdouble c_start = startV;
    tdouble g_1 = LSF_call(c_start,rv_base,rv_prop,betaVec,N_LSF_calls);
  // make sure that root is in selected interval
    tuint c = 0;
    tdouble delta = c_end-c_start;
    bool go_first = true;
    bool go_right = true;
    while ((g_1*g_2>ZERO) && c<iter_max/2) {
      // try to prevent oscillations
        if (go_first) {
          go_right = (g_2<g_1);
          go_first = false;
        } else {
          if (go_right!=(g_2<g_1)) {
            go_first = true;
            delta /= 2;
            const tdouble c_mid = (g_2<g_1)?c_end:c_start;
            c_start = c_mid - delta/2;
            c_end = c_start + delta;
            g_1 = LSF_call(c_start,rv_base,rv_prop,betaVec,N_LSF_calls);
            g_2 = LSF_call(c_end,rv_base,rv_prop,betaVec,N_LSF_calls);
            ++c;
            continue;
          }
        }
      if (g_2<g_1) {
        c_start = c_end;
        c_end += delta;
        g_1 = g_2;
        g_2 = LSF_call(c_end,rv_base,rv_prop,betaVec,N_LSF_calls);
      } else {
        c_end = c_start;
        c_start -= delta;
        g_2 = g_1;
        g_1 = LSF_call(c_start,rv_base,rv_prop,betaVec,N_LSF_calls);
      }
      ++c;
    }
    if (g_1*g_2>ZERO) return (c_end+c_start)/2;
  // start the actual iteration
    tdouble c_res=ZERO, g_3=ZERO;
    while (c<iter_max) {
      c_res = (c_start+c_end)/2;
      g_3 = LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls);
      if ( g_1*g_3>ZERO ) {
        c_start = c_res;
        g_1 = g_3;
      } else {
        c_end = c_res;
        g_2 = g_3;
      }
      if (fabs(g_3)<=tol) {
        found = true;
        break;
      }
      if (fabs(c_end-c_start)<=tol) break;
      ++c;
    }
  // try to deduce information about the topology
    fdright = ( g_1>g_2 );
    if (c_start>c_end) { fdright = !fdright; }  // this should not happend in case of the bisection method
  // maximum number of line-search iterations reached?
    if (c>=iter_max) {
      if (verboseLog) {
        GlobalVar.alert.alert("FlxObjLineSmpl::perform_line_search_bisec","Maximum number of line-search iterations reached.");
      }
    }
  if (fabs(g_3)<=tol) return c_res;
  if (g_1*g_2<=ZERO && c<iter_max) found=true;
  return (c_start+c_end)/2;
}

void FlxObjLineSmpl::hist_push(const tdouble c, const tdouble g)
{
  // search the relevant index
    size_t target_i = line_hist.size()+10;
    for (size_t i=0;i<line_hist.size();++i) {
      if (c<line_hist[i].first) {
        target_i = i;
        break;
      } else if (c==line_hist[i].first) {
        line_hist[i].second = g;
        return;
      }
    }
  // insert pair
    std::pair<tdouble,tdouble> tp(c,g);
    if (target_i>=line_hist.size()) {
      line_hist.push_back(tp);
    } else {
      line_hist.insert(line_hist.begin()+target_i,tp);
    }
}

const tdouble FlxObjLineSmpl::hist_eval(const tdouble betaNorm)
{
  if (line_hist.empty()) return ZERO;
  pdouble res;
  bool is_very_first = true;
  bool is_first = true;
  bool is_failure = false;  // dummy-initialization (to prevent compiler warning)
  tdouble last_c = -tdouble(100);
  for (size_t i=0;i<line_hist.size();++i) {
    tdouble tg = line_hist[i].second;
    if (is_first) {
      is_failure = (tg<=ZERO);
      is_first = false;
    }
    if ( (tg<=ZERO)!=is_failure || tg==ZERO ) {
      const tdouble this_c = (tg==ZERO)?(line_hist[i].first):((line_hist[i].first+line_hist[i-1].first)/2);
      if (is_very_first) {
        if (is_failure) {
          res += rv_Phi(this_c*betaNorm);
        }
        is_very_first = false;
      } else {
        if (is_failure) {
          if (this_c>ZERO && last_c>ZERO) {
            res += rv_Phi(-last_c*betaNorm)-rv_Phi(-this_c*betaNorm);
          } else {
            res += rv_Phi(this_c*betaNorm)-rv_Phi(last_c*betaNorm);
          }
        }
      }
      last_c = this_c;
      if (tg==ZERO) {
        is_first = true;
      } else {
        is_failure = (tg<=ZERO);
        is_first = false;
      }
    }
  };
  // handle ending
    if (is_first) {
      is_failure = (line_hist[line_hist.size()-1].second<=ZERO);
    }
    if (is_failure) {
      if (is_very_first) res += ONE;
      else res += rv_Phi(-last_c*betaNorm);
    }
  return res.cast2double();
}

void FlxObjLineSmpl::task()
{
  tulong Nlsf_calls = 0;
  const tdouble eq_tol = 1e-2;
  // generate the random_creator
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrvsets->eval(true));
    RBRV_constructor RndBox(set_str_vec,data->rbrv_box);
    RndBoxp = &RndBox;
  // evaluate line-search parameters
    const tuint NLS = funNp->cast2tuint();
    const tuint NRV = RndBox.get_NRV();
    const tdouble LS_tol = LS_tolF->cast2positive(true);
    const tuint LS_max_iter = LS_max_iterF->cast2tuint(true);
  // preliminary output
    std::ostream& scout(GlobalVar.slogcout(3));
    scout << "Line sampling: " << std::endl;
    scout << "  number of line-searches:      " << NLS << std::endl;
    scout << "  tolerance of line-search:     " << LS_tol << std::endl;
    scout << "  max. line-search iterations:  " << LS_max_iter << std::endl;
    if (use_bisec) {
    scout << "  line-search algorithm:        " << "bisection method" << std::endl;
    } else {
    scout << "  line-search algorithm:        " << "regula falsi / secant method" << std::endl;
    }
    scout << "  extended line-search:         " << (extended_ls?"yes":"no") << std::endl;
    scout << "  --------------------------------------------------------------------" << std::endl;
  // allocate vectors ...
    flxVec rvy(NRV);
    flxVec rvy_dummy(NRV);
    bool fdright,found;
  // optain the beta-vector
    flxVec betaVec(data->ConstMtxBox.get_Vec(NRV,LS_SPNT->eval(),true),NRV,true);
  // perform an initial line search (to fix the actual tolerance parameter)
    {
      const tdouble c0 = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls,fdright,found,tdouble(0.8),tdouble(1.2));
      scout << "  c0 = " << GlobalVar.Double2String(c0) << std::endl;
      if (found) {
        if (c0>ONE) betaVec *= c0;
      }
    }
    const tdouble betaNorm2 = betaVec.get_Norm2_NOroot();
    const tdouble betaNorm = sqrt(betaNorm2);
  // start with the iteration (line searches)
    // register the progress-bar
      scout << "  preform line searches ...      ";
      FlxProgress prg(scout,!NOTdolog);
      prg.start(NLS);
    const tuint N_LSF_calls_1 = Nlsf_calls;
    pdouble PrEst;
    for (tuint i=0;i<NLS;++i) {
      line_hist.clear();
      pdouble PrEsti;
      // propose a random vector from the prior (standard Normal)
        data->RndCreator.gen_smp(rvy);
      // create in-plane vector
        const tdouble c = -(rvy.operator*(betaVec)/betaNorm2);
        rvy.add(betaVec,c);
      // perform the line search
        const tdouble ci = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls,fdright,found,tdouble(0.8),tdouble(1.2));
        if (found) {
          if (extended_ls) {
            hist_push(ci,ZERO);
          } else {
            if (fdright) {
              PrEsti += rv_Phi(-ci*betaNorm);
            } else {
              PrEsti += rv_Phi(ci*betaNorm);
            }
          }
        }
        if (extended_ls) {
          // search right of 'ci'
            tdouble cinit = ci+2*ONE;
            if (2*ONE>cinit||found==false) cinit = 2*ONE;
            bool fdright2;
            bool found2 = false;
            tdouble ci2 = LSF_call(cinit,rvy,rvy_dummy,betaVec,Nlsf_calls);
            if ( (fdright && ci2>=ZERO) || (!fdright && ci2<=ZERO) ) {
              ci2 = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls,fdright2,found2,ci+(cinit-ci)*0.05,cinit);
              if (found2) {
                hist_push(ci2,ZERO);
              }
            }
          // search left of 'ci'
            if (found) {
              tdouble ci3;
              bool found3=true;
              if (found2 && ci2<ci && fabs(ci2-ci)>eq_tol) {
                ci3 = ci2;
                ci2 *= 2;                // just to make it different from ci3!
              } else {
                cinit = -ONE;
                if (ci-2*ONE<cinit && found) cinit = ci-2*ONE;
                found3=false;
                ci3 = LSF_call(cinit,rvy,rvy_dummy,betaVec,Nlsf_calls);
                if ( (!fdright && ci3>=ZERO) || (fdright && ci3<=ZERO) ) {
                  ci3 = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls,fdright2,found3,cinit,ci-(ci-cinit)*0.05); 
                  if (found3) {
                    hist_push(ci3,ZERO);
                  }
                }
              }
            }
        }
      if (extended_ls) {
        PrEst += hist_eval(betaNorm);
      } else {
        PrEst += PrEsti;
      }
      prg.tick(i+1);
    }
    line_hist.clear();
    // deactivate the progress-bar
      prg.stop();
      scout << "[" << GlobalVar.Double2String(tdouble(Nlsf_calls-N_LSF_calls_1)/tdouble(NLS),false,2) << " calls per line]" << std::endl;
  // determin probability
    theResult = PrEst.cast2double();
    theResult /= NLS;
  // output ...
    GlobalVar.slogcout(4) << " Result of the Integration: " << GlobalVar.Double2String(theResult) << std::endl;
}

FlxObjReadLineSmpl::FlxObjReadLineSmpl()
{
  AllDefParaBox->insert(new FlxOptionalParaFlxString("nataf","sim::rbrvsets",true));
  ParaBox.insert("rbrvsets", "sim::rbrvsets" );
  
  FlxMtxConstFun* mtxvalue = new FlxMtxConstFun( "internal_dummy", NULL );
    AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"ls::spnt"));
    ParaBox.insert("ls_spnt", "ls::spnt" );
  AllDefParaBox->insert(new FlxOptionalParaFun(1e-3,"ls::tol"));
    ParaBox.insert("ls_tol", "ls::tol" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(10),"ls::max_iter"));
    ParaBox.insert("ls_max_iter", "ls::max_iter" );
    
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"ls::extended_ls"));
    ParaBox.insert("extended_ls", "ls::extended_ls" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"ls::use_bisec"));
    ParaBox.insert("use_bisec", "ls::use_bisec" );
}

FlxObjBase* FlxObjReadLineSmpl::read()
{
  reader->getChar('(',false);
  // result-constant
    const std::string cname = reader->getWord(true,false);
    data->ConstantBox.declareC(cname);
    tdouble* theResult = data->ConstantBox.get(cname);
    reader->getChar(';',false);
  // function for number of points
    FlxFunction* funNp = new FlxFunction(funReader,false);
  FlxFunction* fung=NULL;
  try {
      reader->getChar(';',false);
    // function to integrate
      fung = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    read_optionalPara(false); 
    return new FlxObjLineSmpl(get_doLog(), *theResult, funNp, fung, get_optPara_FlxMtxFun("ls_spnt"), get_optPara_FlxFunction("ls_tol"), get_optPara_FlxFunction("ls_max_iter"), get_optPara_bool("extended_ls"), get_verboseLog(), get_optPara_FlxString("rbrvsets"), get_optPara_bool("use_bisec") );
  } catch (FlxException& e) {
    if (funNp) delete funNp;
    if (fung) delete fung;
    throw;
  }
}



SuS_csm_evalStorage::SuS_csm_evalStorage(FlxFunction* h_value, FlxString* kernel_name, FlxString* MCMCmethS, FlxFunction* csm_p, FlxFunction* csm_nmax, FlxFunction* csm_p_single, FlxFunction* csm_nmax_single, FlxFunction* dcs_pSD)
: h_value(h_value), kernel_name(kernel_name), MCMCmethS(MCMCmethS), 
  csm_p(csm_p), csm_nmax(csm_nmax), csm_p_single(csm_p_single), csm_nmax_single(csm_nmax_single),
  dcs_pSD(dcs_pSD)
{
  
}

SuS_csm_evalStorage::~SuS_csm_evalStorage()
{
  delete h_value;
  delete kernel_name;
  delete MCMCmethS;
  delete csm_p;
  delete csm_nmax;
  delete csm_p_single;
  delete csm_nmax_single;
  delete dcs_pSD;
}

FlxBayUP_csm_base* SuS_csm_evalStorage::eval(FlxBayUp_Update_List* list)
{
  FlxFunction* h_copy = NULL;
  FlxBayUP_csm_base* csm = NULL;
  try {
    // get MCMC meth
      std::string mStr = MCMCmethS->eval_word(true);
    const tdouble h = h_value->cast2positive(false);
    if (!(list->get_adpt_ctrl().is_adaptive())) {
      if (h_value->dependOn_Const(data->ConstantBox.get("sus_iter"))) {
        h_copy = new FlxFunction(*h_value);
      }
    }
    FlxRndCreator& RndCreator = list->parent.updater.get_RndCreator();
    if (mStr=="cwmh") {
      mStr = kernel_name->eval_word(true);
      csm = new FlxBayUP_csm_cwmh_MCMC(RndCreator,mStr,h,h_copy);
    } else if (mStr=="cov") {
      mStr = kernel_name->eval_word(true);
      const tdouble csm_pv = csm_p->cast2positive(false);
      const tuint csm_nmaxv = csm_nmax->cast2tuintW0(false);
      const tdouble csm_p_singlev = csm_p_single->cast2positive(false);
      const tuint csm_nmax_singlev = csm_nmax_single->cast2tuintW0(false);
      csm = new FlxBayUP_csm_cov_MCMC(RndCreator,list->get_Nrv(),mStr,h,h_copy,csm_pv,csm_nmaxv,csm_p_singlev,csm_nmax_singlev,*list);
    } else if (mStr=="csusmh") {
      csm = new FlxBayUP_csm_csus_MCMC(RndCreator,h,h_copy);
    } else if (mStr=="dcs") {
      const tdouble dcs_pSDv = dcs_pSD->cast2positive_or0(false);
      csm = new FlxBayUP_csm_dcs_MCMC(RndCreator,h,dcs_pSDv,h_copy,*list);
    } else if (mStr=="tmcmc") {
      csm = new FlxBayUP_csm_TMCMC(RndCreator,list->get_Nrv(),h,h_copy);
    } else {
      std::ostringstream ssV;
      ssV << "Unknown ID for an MCMC method (" << mStr << ").";
      throw FlxException("SuS_csm_evalStorage::eval",ssV.str());
    }
    csm->register_adpt_ctrl(&(list->get_adpt_ctrl()));
  } catch (FlxException &e) {
    delete list;
    if (h_copy) delete h_copy;
    throw;
  }
  return csm;
}

FlxObjSuS::FlxObjSuS(const bool dolog, const std::string& ostreamV, FlxFunction* Nc, FlxFunction* Ncl, FlxFunction* max_runs, const FlxBayUp_Update_List::randomizeTec randomize, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const Flx_SuS_Control& susControl, SuS_csm_evalStorage* csm_eval, FlxString* rbrvsets, FlxFunction* lsf)
: FlxObjOutputBase(dolog,ostreamV),Nc(Nc),Ncl(Ncl),max_runs(max_runs), randomize(randomize), adpt_ctrl(adpt_ctrl), 
  susControl(susControl), csm_eval(csm_eval), rbrvsets(rbrvsets), lsf(lsf)
{

}

FlxObjSuS::~FlxObjSuS()
{
  delete Nc;
  delete Ncl;
  delete max_runs;
  delete adpt_ctrl;
  delete csm_eval;
  if (rbrvsets) delete rbrvsets;
  if (lsf) delete lsf;
}

void FlxObjSuS::task()
{
  if (rbrvsets==NULL || lsf==NULL) throw FlxException_Crude("FlxObjSuS::task");
  const std::string setstr = rbrvsets->eval(true);
  if (FlxObjReadSuS::lastSuS) delete FlxObjReadSuS::lastSuS;
  FlxObjReadSuS::lastSuS = new flxBayUp(setstr);
  // define LSF-object
    FlxObjReadSuS::lastSuS->set_globalLkl(*lsf,false,flxBayUp::RA);
  const tuint NcV = Nc->cast2tuint(false);
  const tuint NclV = Ncl->cast2tuint(false);
  const tuint mr = max_runs->cast2tuint(false);
  const tuint NsfV = NcV*NclV;
  FlxBayUp_Update_List* list = new FlxBayUp_Update_List(*(FlxObjReadSuS::lastSuS),NcV,NclV,NsfV,0,randomize,adpt_ctrl->copy(),mr,false,FlxBayUp_Update_List::RASUBSIM,false,susControl.find_multiples);
  FlxBayUP_csm_base* csm = csm_eval->eval(list);
  FlxObjReadSuS::lastSuS->updater.update(list,csm,false,susControl);
}

FlxObjReadSuS::FlxObjReadSuS()
{
  FlxBayUp_Update::define_constants();
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(50.),"sus::max_runs"));
    ParaBox.insert("max_runs", "sus::max_runs" );
  AllDefParaBox->insert(new FlxOptionalParaStream("rndpick","sus::randomize"));
    ParaBox.insert("randomize", "sus::randomize" );
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"sus::prt_alert"));
    ParaBox.insert("prt_alert", "sus::prt_alert" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"sus::ext_out"));
    ParaBox.insert("ext_out", "sus::ext_out" );
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"sus::comp_gamma"));
    ParaBox.insert("comp_gamma", "sus::comp_gamma" );
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"sus::consider_seed_corr"));
    ParaBox.insert("consider_seed_corr", "sus::consider_seed_corr" );
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"sus::consider_pi_corr"));
    ParaBox.insert("consider_pi_corr", "sus::consider_pi_corr" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"sus::empirical_corr"));
    ParaBox.insert("empirical_corr", "sus::empirical_corr" );
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"sus::find_multiples"));
    ParaBox.insert("find_multiples", "sus::find_multiples" );
  AllDefParaBox->insert(new FlxOptionalParaFlxString("","sus::write_smpls",false));
    ParaBox.insert("write_smpls", "sus::write_smpls" );
  // ------------------------------------------------------
  // BUS
  // ------------------------------------------------------
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"bus::pa_maxl"));
    ParaBox.insert("pa_maxl", "bus::pa_maxl" );
  // ------------------------------------------------------
  // TMCMC
  // ------------------------------------------------------
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"tmcmc::target_cov"));
    ParaBox.insert("target_cov", "tmcmc::target_cov" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"tmcmc::tmcmc_update_weights"));
    ParaBox.insert("tmcmc_update_weights", "tmcmc::tmcmc_update_weights" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"tmcmc::tmcmc_alpha"));
    ParaBox.insert("tmcmc_alpha", "tmcmc::tmcmc_alpha" );
  // ------------------------------------------------------
  // Line sampling
  // ------------------------------------------------------
    FlxMtxConstFun* mtxvalue = new FlxMtxConstFun( "internal_dummy", NULL );
    AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"ls::spnt"));
    ParaBox.insert("ls_spnt", "ls::spnt" );
  AllDefParaBox->insert(new FlxOptionalParaFun(1e-3,"ls::tol"));
    ParaBox.insert("ls_tol", "ls::tol" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(10),"ls::max_iter"));
    ParaBox.insert("ls_max_iter", "ls::max_iter" );
  // ------------------------------------------------------
  // credible intervals
  // ------------------------------------------------------
    AllDefParaBox->insert(new FlxOptionalParaStream("none","sus::credest"));
      ParaBox.insert("credest", "sus::credest" );
      std::vector<FlxFunction*> vecV;
      vecV.push_back(new FlxFunction(new FunNumber(0.90)));
      vecV.push_back(new FlxFunction(new FunNumber(0.95)));
      FlxObjBase* blockV = new FlxObjMtxConstNewU(false, new FlxMtxConstFun("internal_susrelpc"), vecV, 2, 1);
      mtxvalue = new FlxMtxConstFun( "internal_susrelpc", blockV );
      AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"sus::credible"));
      ParaBox.insert("credible", "sus::credible" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(1e5),"sus::n_smpl_cred"));
      ParaBox.insert("n_smpl_cred", "sus::n_smpl_cred" );
  // ------------------------------------------------------
  // MCMC methods
  // ------------------------------------------------------
  AllDefParaBox->insert(new FlxOptionalParaFlxString("cwmh","sus::mcmc_algo",true));
    ParaBox.insert("mcmc_algo", "sus::mcmc_algo" );
  AllDefParaBox->insert(new FlxOptionalParaFlxString("gauss","sus::kernel",true));
    ParaBox.insert("kernel", "sus::kernel" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"sus::kernel_h"));
    ParaBox.insert("kernel_h", "sus::kernel_h" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.2),"sus::csm_p"));
    ParaBox.insert("csm_p", "sus::csm_p" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(5e3),"sus::csm_nmax"));
    ParaBox.insert("csm_nmax", "sus::csm_nmax" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.2),"sus::csm_p_single"));
    ParaBox.insert("csm_p_single", "sus::csm_p_single" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(10.),"sus::csm_nmax_single"));
    ParaBox.insert("csm_nmax_single", "sus::csm_nmax_single" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(ZERO),"sus::dcs_psd"));
    ParaBox.insert("dcs_psd", "sus::dcs_psd" );
  // ------------------------------------------------------
  // Adaptive settings
  // ------------------------------------------------------
  AllDefParaBox->insert(new FlxOptionalParaFlxString("linear","sus::adaptive_meth",true));
    ParaBox.insert("adaptive_meth", "sus::adaptive_meth" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(90.),"sus::adaptive_afternmodruns"));
    ParaBox.insert("adaptive_afternmodruns", "sus::adaptive_afternmodruns" );
  AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(2.),"sus::adaptive_smpl_order"));
    ParaBox.insert("adaptive_smpl_order", "sus::adaptive_smpl_order" );
  // adaptive_meth = linear
    AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"sus::adaptive_factor"));
      ParaBox.insert("adaptive_factor", "sus::adaptive_factor" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.3),"sus::adaptive_lower"));
      ParaBox.insert("adaptive_lower", "sus::adaptive_lower" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.5),"sus::adaptive_upper"));
      ParaBox.insert("adaptive_upper", "sus::adaptive_upper" );
  // adaptive_meth = log
    AllDefParaBox->insert(new FlxOptionalParaFun(data->ReadManager.parse_function("1/sqrt(sus_iadpt)"),"sus::adaptive_f1"));
      ParaBox.insert("adaptive_f1", "sus::adaptive_f1" );
    AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"sus::adaptive_f2"));
      ParaBox.insert("adaptive_f2", "sus::adaptive_f2" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.44),"sus::adaptive_acr"));
      ParaBox.insert("adaptive_acr", "sus::adaptive_acr" );
  // adaptive_meth = dcs
    AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"sus::adaptive_dcsf"));
      ParaBox.insert("adaptive_dcsf", "sus::adaptive_dcsf" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.5),"sus::psd_max"));
      ParaBox.insert("psd_max", "sus::psd_max" );
  // adaptive_meth = velo
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.5),"sus::adaptive_vspread"));
      ParaBox.insert("adaptive_vspread", "sus::adaptive_vspread" );
  // adaptive_meth = opti_jump
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.2),"sus::adaptive_acr_min"));
      ParaBox.insert("adaptive_acr_min", "sus::adaptive_acr_min" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.95),"sus::adaptive_esjd_scale"));
      ParaBox.insert("adaptive_esjd_scale", "sus::adaptive_esjd_scale" );
    AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(1e4),"sus::adaptive_nmax"));
      ParaBox.insert("adaptive_nmax", "sus::adaptive_nmax" );
    AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"sus::adaptive_pw_p1"));
      ParaBox.insert("adaptive_pw_p1", "sus::adaptive_pw_p1" );
    AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"sus::adaptive_pw_p2"));
      ParaBox.insert("adaptive_pw_p2", "sus::adaptive_pw_p2" );
    AllDefParaBox->insert(new FlxOptionalParaFun(0.01,"sus::adaptive_aeps"));
      ParaBox.insert("adaptive_aeps", "sus::adaptive_aeps" );
  // only for Subset simulation - not for Baysesian updating
    AllDefParaBox->insert(new FlxOptionalParaFlxString("nataf","sim::rbrvsets",true));
      ParaBox.insert("rbrvsets", "sim::rbrvsets" );
}

FlxObjReadSuS::~FlxObjReadSuS()
{
  if (lastSuS) {
    delete lastSuS;
    lastSuS = NULL;
  }
}

FlxBayUp_Update_List::randomizeTec FlxObjReadSuS::get_randomize_id()
{
  const std::string rndTecStr = get_optPara_string("randomize",true);
  if (rndTecStr=="none") {
    return FlxBayUp_Update_List::NONE;
  } else if (rndTecStr=="init") {
    return FlxBayUp_Update_List::INIT;
  } else if (rndTecStr=="rndpick") {
    return FlxBayUp_Update_List::RNDPICK;
  } else {
    std::ostringstream ssV;
    ssV << "Unknown randomization technique (" << rndTecStr << ").";
    throw FlxException_NeglectInInteractive("FlxObjReadBayUp_update::read_3", ssV.str() );
  }
}

flxBayUp_adaptive_ctrl_base* FlxObjReadSuS::get_adpt_ctrl()
{
  // smpl_order
    const tuint smpl_order = get_optPara_tuint_from_FlxFunction("adaptive_smpl_order",true,false);
  // adaptive scheme
    std::string asstr = get_optPara_word_from_FlxString("adaptive_meth",true);
    if (asstr=="linear") {
      return new flxBayUp_adaptive_ctrl_bounds(get_optPara_FlxFunction("adaptive_factor"),get_optPara_FlxFunction("adaptive_lower"),get_optPara_FlxFunction("adaptive_upper"),
              get_optPara_FlxFunction("adaptive_afternmodruns"), smpl_order);
    } else if (asstr=="log") {
      return new flxBayUp_adaptive_ctrl_log(get_optPara_FlxFunction("adaptive_f1"),get_optPara_FlxFunction("adaptive_f2"),get_optPara_FlxFunction("adaptive_acr"),get_optPara_FlxFunction("adaptive_afternmodruns"), smpl_order);
    } else if (asstr=="dcs") {
      return new flxBayUp_adaptive_ctrl_dcs(get_optPara_FlxFunction("adaptive_afternmodruns"), get_optPara_FlxFunction("adaptive_dcsf"), get_optPara_FlxFunction("psd_max"), smpl_order);
    } else if (asstr=="velo") {
      return new flxBayUp_adaptive_ctrl_velo(data->RndCreator, get_optPara_FlxFunction("adaptive_vspread"), get_optPara_FlxFunction("adaptive_afternmodruns"), smpl_order);
    } else if (asstr=="opti_jump") {
      return new flxBayUp_adaptive_ctrl_opti_jump(data->RndCreator, get_optPara_FlxFunction("adaptive_acr_min"), get_optPara_FlxFunction("adaptive_esjd_scale"), get_optPara_FlxFunction("adaptive_pw_p1"), get_optPara_FlxFunction("adaptive_pw_p2"), get_optPara_FlxFunction("adaptive_aeps"), get_optPara_FlxFunction("adaptive_Nmax"), get_optPara_FlxFunction("adaptive_afternmodruns"), smpl_order);
    } else {
      std::ostringstream ssV;
      ssV << "Unknown ID (" << asstr << ") for adaptive_meth.";
      throw FlxException_NeglectInInteractive("FlxObjReadSuS::get_adpt_ctrl_3", ssV.str() );
    }
}

SuS_csm_evalStorage* FlxObjReadSuS::get_csm_eval()
{
  return new SuS_csm_evalStorage(
    get_optPara_FlxFunction("kernel_h"),get_optPara_FlxString("kernel"),get_optPara_FlxString("mcmc_algo"),
    get_optPara_FlxFunction("csm_p"),get_optPara_FlxFunction("csm_nmax"),get_optPara_FlxFunction("csm_p_single"),get_optPara_FlxFunction("csm_nmax_single"),
    get_optPara_FlxFunction("dcs_psd")
  );
}

const Flx_SuS_Control FlxObjReadSuS::get_susControl()
{
  Flx_SuS_Control susControl;
    susControl.prt_alert = get_optPara_bool("prt_alert");
    susControl.verbose = get_optPara_bool("ext_out");
    susControl.credEst = Flx_SuS_Control::parse_credibleEstim(get_optPara_string("credest",true));
    susControl.pc = get_optPara_FlxMtxFun("credible");
    susControl.consider_seed_corr = get_optPara_bool("consider_seed_corr");
    susControl.consider_pi_corr = get_optPara_bool("consider_pi_corr");
      if (!susControl.consider_seed_corr) susControl.consider_pi_corr = false;
    susControl.empirical_corr = get_optPara_bool("empirical_corr");
    susControl.N_cred_smpl = get_optPara_tuint_from_FlxFunction("n_smpl_cred",false);
    susControl.find_multiples = get_optPara_bool("find_multiples");
    susControl.os_samples = get_optPara_FlxString("write_smpls");
    susControl.TMCMC_target_COV = get_optPara_FlxFunction("target_cov");
    susControl.TMCMC_update_weights = get_optPara_tuint_from_FlxFunction("tmcmc_update_weights",true);
    susControl.TMCMC_alpha = get_optPara_FlxFunction("tmcmc_alpha");
    susControl.LS_SPNT = get_optPara_FlxMtxFun("ls_spnt");
    susControl.LS_tol = get_optPara_FlxFunction("ls_tol");
    susControl.LS_max_iter = get_optPara_FlxFunction("ls_max_iter");
    susControl.pa_maxL = get_optPara_FlxFunction("pa_maxl");
  if (susControl.verbose || (susControl.credEst!=Flx_SuS_Control::none && susControl.credEst!=Flx_SuS_Control::simpleBayes) ) {
    susControl.comp_gamma = true;
  }
  return susControl;
}

FlxObjBase* FlxObjReadSuS::read () {
  FlxFunction* funNcl=NULL; FlxFunction* funNc=NULL; FlxFunction* fung=NULL;
  flxBayUp_adaptive_ctrl_base* adpt_ctrl = NULL;
  SuS_csm_evalStorage* csm_eval = NULL;
  try {
    reader->getChar('(',false);
    // number of chains
      funNc = new FlxFunction(funReader,false);
      reader->getChar(';',false);
    // number of samples
      funNcl = new FlxFunction(funReader,false);
      reader->getChar(';',false);
    // limit-state function
      fung = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    
    FlxBayUp_Update_List::randomizeTec rnd_id = get_randomize_id();
    adpt_ctrl = get_adpt_ctrl();
    csm_eval = get_csm_eval();
    
    Flx_SuS_Control susControl = get_susControl();
    
    return new FlxObjSuS(get_doLog(),get_stream(),funNc,funNcl,get_optPara_FlxFunction("max_runs"),
                                  rnd_id,adpt_ctrl,susControl,csm_eval,get_optPara_FlxString("rbrvsets"),fung);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSubSetSim::read_99",1);
    if (funNcl) delete funNcl;
    if (funNc) delete funNc;
    if (fung) delete fung;
    if (adpt_ctrl) delete adpt_ctrl;
    if (csm_eval) delete csm_eval;
    throw;
  }
}

void FlxObjSus_level_info::task()
{
  const tuint pid = pidf->cast2tuint(false);
  const tuint pid2 = pidf2?(pidf2->cast2tuintW0(false)):0;
  const std::string vecs = VecStr->eval();
  if (nameID) {
    const std::string name = nameID->eval_word(true);
    flxBayUp& bu = BayUpBox->get(name);
    bu.updater.get_sus_level_info(vecs,pid,pid2);
  } else {
    if (FlxObjReadSuS::lastSuS==NULL) {
      throw FlxException("FlxObjSus_level_info::task","You must execute Subset Simulation before you can call 'sus_level_info'.");
    }
    FlxObjReadSuS::lastSuS->updater.get_sus_level_info(vecs,pid,pid2);
  }
}

FlxObjBase* FlxObjReadSus_level_info::read()
{
  FlxMtxConstFun* VecStr = new FlxMtxConstFun(false);
  FlxString* nameID = NULL;
  FlxFunction* pidf = NULL;  
  FlxFunction* pid2 = NULL;  
  try {
    reader->getChar('=',false);
    if (reader->whatIsNextChar()==':') {
      reader->getExpr("::",false);
    } else {
      nameID = new FlxString(false,false);
    }
    reader->getChar('(',false);
    pidf = new FlxFunction(funReader,false);
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      pid2 = new FlxFunction(funReader,false);
    }
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjSus_level_info(get_doLog(),VecStr,nameID,pidf,pid2);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSus_level_info::read",1);
    delete VecStr;
    if (nameID) delete nameID;
    if (pidf) delete pidf;
    if (pid2) delete pid2;
    throw;
  }
}

void FlxObjFORM_betaSensitivities::task()
{
  const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rvsets->eval(true));
  RBRV_constructor RndBox(set_str_vec,data->rbrv_box);
  tuint xN = RndBox.get_NRV();
  flxVec rvyp(data->ConstMtxBox.get_Vec(rvy->eval(),xN),xN,true);
  flxVec svp(data->ConstMtxBox.get_Vec(sv->eval(),xN),xN);
  RndBox.set_smp(rvyp);
  FlxObjFORM::sensitivities(rvyp,RndBox,sout(),&svp);
}

FlxObjBase* FlxObjReadFORMbetaSensitivities::read()
{
  FlxMtxConstFun* sv = new FlxMtxConstFun(false);
  FlxMtxConstFun* rvy = NULL;
  FlxString* rvsets = NULL;
  try {
    reader->getChar('=',false);
    rvy = new FlxMtxConstFun(true);
    reader->getChar('(',false);
    rvsets = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjFORM_betaSensitivities(get_doLog(),get_stream(),sv,rvy,rvsets);
  } catch (FlxException &e) {
    delete sv;
    if (rvy) delete rvy;
    if (rvsets) delete rvsets;
    throw;
  }
}

FlxObjReadFORM_base::FlxObjReadFORM_base()
{
  // optional parameters
    // maxiter
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(100),"form::maxiter"));
      ParaBox.insert("maxiter", "form::maxiter" );
    // fdstep
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(1e-6),"form::fdstep"));
      ParaBox.insert("fdstep", "form::fdstep" );
    // fdstep
      AllDefParaBox->insert(new FlxOptionalParaFun(2*ONE,"form::epsdyf"));
      ParaBox.insert("epsdyf", "form::epsdyf" );
    // eps1
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(1e-5),"form::eps1"));
      ParaBox.insert("eps1", "form::eps1" ); 
    // eps2
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(1e-5),"form::eps2"));
      ParaBox.insert("eps2", "form::eps2" );
    // ihlrf_lambda_start
      AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"form::ihlrf_lambda_start"));
      ParaBox.insert("ihlrf_lambda_start", "form::ihlrf_lambda_start" );
    // ihlrf_epsilon
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(0.4),"form::ihlrf_epsilon"));
      ParaBox.insert("ihlrf_epsilon", "form::ihlrf_epsilon" );
    // ihlrf_reduce
      AllDefParaBox->insert(new FlxOptionalParaFun(ONE/2,"form::ihlrf_reduce"));
      ParaBox.insert("ihlrf_reduce", "form::ihlrf_reduce" );
    // xstart
      FlxMtxConstFun* mtxvalue = new FlxMtxConstFun( "internal_formxstart" );
      AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"form::xstart"));
      ParaBox.insert("xstart", "form::xstart" );
    // dx_min
      mtxvalue = new FlxMtxConstFun( "internal_formdxmin" );
      AllDefParaBox->insert(new FlxOptionalParaMtxFun(mtxvalue,"form::dxmin"));
      ParaBox.insert("dxmin", "form::dxmin" ); 
    // dxdyanalytical
      AllDefParaBox->insert(new FlxOptionalParaBool(true,"form::dxdyanalytical"));
      ParaBox.insert("dxdyanalytical", "form::dxdyanalytical" );
    // fd_method
      AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"form::fd_method"));
      ParaBox.insert("fd_method", "form::fd_method" );
    // opt_method
      AllDefParaBox->insert(new FlxOptionalParaFun(2*ONE,"form::opt_method"));
      ParaBox.insert("opt_method", "form::opt_method" );
}

FlxObjReadFORM_pdf::FlxObjReadFORM_pdf()
{
  // optional parameters
    // lbound
      AllDefParaBox->insert(new FlxOptionalParaFun(tdouble(100),"form_pdf::intervals"));
      ParaBox.insert("intervals", "form_pdf::intervals" );
    // Default verbose-setting
      AllDefParaBox->insert(new FlxOptionalParaBool(false,"form_pdf::verbose"));
      ParaBox.insert("vlog", "form_pdf::verbose" );
}

FlxObjBase* FlxObjReadFORM_pdf::read()
{
  reader->getChar('(',false);
  FlxFunction* rvfun = new FlxFunction(funReader,false);
  FlxFunction *lbound=NULL, *ubound=NULL;
  try {
    reader->getChar(';',false);
    reader->getChar('[',false);
    lbound = new FlxFunction(funReader,false);
    reader->getChar(';',false);
    ubound = new FlxFunction(funReader,false);
    reader->getChar(']',false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjFORM_pdf(get_doLog(),rvfun,lbound,ubound,get_optPara_FlxFunction("intervals"),get_optPara_FlxMtxFun("xstart"),get_optPara_FlxFunction("fdstep"), get_optPara_FlxFunction("epsdyf"), get_optPara_FlxFunction("eps1"), get_optPara_FlxFunction("eps2"), get_optPara_FlxFunction("ihlrf_lambda_start"), get_optPara_FlxFunction("ihlrf_epsilon"), get_optPara_FlxFunction("ihlrf_reduce"), get_optPara_tuint_from_FlxFunction("maxiter",false),get_optPara_bool("vlog"),get_stream(),get_optPara_bool("dxdyanalytical"), get_optPara_FlxMtxFun("dxmin"), get_optPara_int_from_FlxFunction("fd_method"), get_optPara_int_from_FlxFunction("opt_method"), get_optPara_FlxString("rbrvsets") );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFORM_pdf::read_1",1);
    delete rvfun;
    if (lbound) delete lbound;
    if (ubound) delete ubound;
    throw;
  }
}

FlxObjReadFORM::FlxObjReadFORM(const bool only_partial) : only_partial(only_partial)
{
  // optional parameters
    // betaDP -> name of const variable where to store
      AllDefParaBox->insert(new FlxOptionalParaFlxString("","form::betadp",false));
      ParaBox.insert("betadp", "form::betadp" ); 
    // Default verbose-setting
      AllDefParaBox->insert(new FlxOptionalParaBool(true,"flxlog::verbose"));
      ParaBox.insert("vlog", "flxlog::verbose" );
}

FlxObjBase* FlxObjReadFORM::read()
{
  reader->getChar('(',false);
  // result-constant
    const std::string cnamey = reader->getWord(true,false);
    data->ConstMtxBox.declareC(cnamey);
    std::string cnamex;
    if (!only_partial) {
      reader->getChar(';',false);
      cnamex = reader->getWord(true,false);
      data->ConstMtxBox.declareC(cnamex);
    }
    reader->getChar(';',false);
  FlxFunction* lsf = new FlxFunction(funReader,false);
  FlxString* rbrvsets = NULL;
  try {
    reader->getChar(';',false);
    rbrvsets = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjFORM( get_doLog(), cnamey, cnamex, get_optPara_FlxMtxFun("xstart"), lsf, get_optPara_FlxFunction("fdstep"), get_optPara_FlxFunction("epsdyf"), get_optPara_FlxFunction("eps1"), get_optPara_FlxFunction("eps2"), get_optPara_FlxFunction("ihlrf_lambda_start"), get_optPara_FlxFunction("ihlrf_epsilon"), get_optPara_FlxFunction("ihlrf_reduce"), get_optPara_tuint_from_FlxFunction("maxiter",false), get_optPara_FlxString("betadp"), get_optPara_bool("vlog"), get_optPara_bool("dxdyanalytical"), get_optPara_FlxMtxFun("dxmin"), get_optPara_int_from_FlxFunction("fd_method"), get_optPara_int_from_FlxFunction("opt_method"), rbrvsets, only_partial );
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadFORM::read",1);
    delete lsf;
    if (rbrvsets) delete rbrvsets;
    throw;
  }
}

FlxObjFORM_base::~FlxObjFORM_base()
{
  if (LSF) delete LSF;
  delete fdstep;
  delete epsdyfF;
  delete eps1;
  delete eps2;
  delete iHLRF_lambda_start;
  delete iHLRF_epsilon;
  delete iHLRF_reduce;
  delete xstart;
  delete dx_min;
  delete rbrvsets;
  if (RndBox) delete RndBox;
}

void FlxObjFORM_base::update_Start()
{
  // generate the random_creator
    if (RndBox) delete RndBox;
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(rbrvsets->eval(true));
    RndBox = new RBRV_constructor(set_str_vec,data->rbrv_box);
  // get DIM
    if (RndBox->get_NRV()!=RndBox->get_NOX()) {
      std::ostringstream ssV;
      ssV << "Number of random variables in original space does not equal number of random variables in standard normal space.";
      throw FlxException("FlxObjFORM_base::update_Start_2", ssV.str() );
    }
    DIM = RndBox->get_NRV();
    if (DIM==0) {
      std::ostringstream ssV;
      ssV << "FORM cannot be executed because no random variables exist.";
      throw FlxException("FlxObjFORM_base::update_Start_2", ssV.str() );
    }
}

void FlxObjFORM_base::eval_xStart(flxVec& x)
{
  tuint xN = 0;
  const std::string xsn = xstart->eval();
  if (xsn=="internal_formxstart") {
    RndBox->get_mean_Vec(x.get_tmp_vptr());
  } else {  
    flxVec xp(data->ConstMtxBox.get_Vec(xsn,xN),xN);
    if (xN != DIM) {
      std::ostringstream ssV;
      ssV << "Vector sizes do not match.\n\tDIM_required=" << DIM << "; DIM_x=" << xN;
      throw FlxException("FlxObjFORM_base::eval_xStart", ssV.str() );
    }
    x = xp;
  }
}

flxVec FlxObjFORM_base::do_FORM(flxVec& x, flxVec& y, tdouble& beta_new, tuint& LSFcalls, const bool only_partial)
{
  const tdouble sqrtEps = GlobalVar.sqrtEps;        // square root of machine precision
  // check parameters
    GlobalVar.slogcout(4) << "  fd_method:          " << fd_method << " - ";
    switch (fd_method) {
      case 1:
        GlobalVar.slogcout(4) << "forward difference";
        break;
      case 2:
        GlobalVar.slogcout(4) << "central difference";
        break;
      case 3:
        GlobalVar.slogcout(4) << "backward difference";
        break;
      default:
        std::ostringstream ssV;
        ssV << "Finite difference method with ID '" << fd_method << "' not known.";
        throw FlxException("FlxObjFORM_base::do_FORM_1", ssV.str() );
    }
    if (!only_partial) {
      GlobalVar.slogcout(4) << std::endl << "  opt_method:         " << opt_method << " - ";
      switch (opt_method) {
        case 1:
          GlobalVar.slogcout(4) << "HLRF-method";
          break;
        case 2:
          GlobalVar.slogcout(4) << "iHLRF-method";
          break;
        case 3:
          throw FlxException_NotImplemented("FlxObjFORM_base::do_FORM_2");
          break;
        default:
          std::ostringstream ssV;
          ssV << "Optimization algorithm with ID '" << fd_method << "' not known.";
          throw FlxException("FlxObjFORM_base::do_FORM_3", ssV.str() );
      }
    }
    const tdouble epsdyf = epsdyfF->cast2positive();
    if (epsdyf<ONE) {
      GlobalVar.alert.alert("FlxObjFORM_base::do_FORM_4","'epsdyf' should not be smaller than 1.0.");
    }
    GlobalVar.slogcout(4) << std::endl;
    GlobalVar.slogcout(4) <<   "  fdstep:             " << fdstep->write()  << "\t(" << GlobalVar.Double2String(sqrtEps) << ")" << std::endl;
    GlobalVar.slogcout(4) <<   "  epsdyf:             " << GlobalVar.Double2String(epsdyf) << std::endl;
    GlobalVar.slogcout(4) <<   "  dxmin:              " << dx_min->write() << std::endl;
    if (!only_partial) {
      GlobalVar.slogcout(4) <<   "  eps1:               " << eps1->write() << std::endl;
      GlobalVar.slogcout(4) <<   "  eps2:               " << eps2->write() << std::endl;
      GlobalVar.slogcout(4) <<   "  maxiter:            " << maxIter << std::endl;
      if (opt_method==2) {
        GlobalVar.slogcout(4) << "  iHLRF_epsilon:      " << iHLRF_epsilon->write() << std::endl;
        GlobalVar.slogcout(4) << "  iHLRF_lambda_start: " << iHLRF_lambda_start->write() << std::endl;
        GlobalVar.slogcout(4) << "  iHLRF_reduce:       " << iHLRF_reduce->write() << std::endl;
      }
    }
  // evaluate dxmin
    flxVec dxmin(DIM);
    {
      const std::string dxmns = dx_min->eval();
      if (dxmns=="internal_formdxmin") {
        RndBox->get_sd_Vec(dxmin.get_tmp_vptr());
        dxmin *= tdouble(1e-6);
      } else {
        const flxVec dxminhlp(data->ConstMtxBox.get_Vec(DIM,dxmns,true),DIM);
        dxmin = dxminhlp;
      }
    }
  // calculate z_0  (value of LSF at y=0)
    tdouble z_0, z_old;
    flxVec d(y.get_N());                // iHLRF direction vector
    flxVec y_new(y.get_N());        // iHLRF new position vector
    flxVec xhelp(DIM);
    d.set_zero();
    RndBox->set_smp(d);
    z_0 = LSF->calc(); ++LSFcalls;
    const bool isNeg_z_0 = (z_0<ZERO);
    GlobalVar.slogcout(4) << "  lsf(y=0)=" << format("%9.2e") % z_0 << std::endl;
    // compute start-point
      RndBox->get_x_Vec(xhelp.get_tmp_vptr());
      xhelp -= x;
      if (xhelp.get_NormMax()<=GlobalVar.TOL()) {
        z_old = z_0;
        y = d;
      } else {
        RndBox->set_smp_x_transform(x);
        z_old = LSF->calc(); ++LSFcalls;
        RndBox->get_y_Vec(y.get_tmp_vptr());
      }
    // make sure z_0 is large enough
      if (isNeg_z_0) z_0 = fabs(z_0);
      if (z_0 < epsdyf*sqrtEps) {
        z_0 = epsdyf*sqrtEps;
      }
    z_old /= z_0;
  // Write log-message
    if (verboseLog) {
      GlobalVar.slog(4) << " Start point: x=" << x << "; y=" << y << ";" << std::endl;
    }
    tuint loopc = 0;
  beta_new = y.get_Norm2();
  
  if (!only_partial) {
    GlobalVar.slogcout(4) << "  Iter  0: lsf_err=" << format("%9.2e") % z_old << std::endl;
  }
  flxVec dzdy(DIM);
  tdouble beta_old, zo, zu,fd,fd_max,dx,lambda=-9999.;
  FlxMtxLTri dxdy(DIM);
  do {
    beta_old = beta_new;
    fd = fdstep->calc();
    if (fd<sqrtEps) {
      fd = sqrtEps;
      GlobalVar.alert.alert("FlxObjFORM_base::do_FORM_5","Set finite differenc step to machine precision.");
    }
    fd_max = fd;
    
    if (dxdyAnalytical) {
      // calc dx/dy
        RndBox->calc_Jinv(dxdy);
      // calculate dz/dx
      for (tuint i = 0; i < DIM; ++i) {
        xhelp = x;
        // choose an appropriate step-size for finite-difference
          // check y-value (machine precision)
            if (fabs((fd*x[i])) < fabs(epsdyf*sqrtEps*z_old)) {        // z determins step-size
              if (fd_method!=3) {        // forward (central)
                xhelp[i] += fabs(epsdyf*sqrtEps*z_old);
                dx = xhelp[i]-x[i];
              } else {                        // backwards
                xhelp[i] -= fabs(epsdyf*sqrtEps*z_old);
                dx = x[i]-xhelp[i];
              }
              const tdouble tr = dx/x[i];
              if (tr>fd_max) fd_max = tr;
            } else {                                                        // x determins step-size
              if (fd_method!=3) {        // forward (central)
                xhelp[i] = (x[i]>=ZERO)?((ONE+fd)*x[i]):((ONE-fd)*x[i]);
                dx = xhelp[i]-x[i];
              } else {                        // backwards
                xhelp[i] = (x[i]>=ZERO)?((ONE-fd)*x[i]):((ONE+fd)*x[i]);
                dx = x[i]-xhelp[i];
              }
            }
          // ensure dx>dxmin
            if (dx<dxmin[i]) {
              if (fd_method!=3) {        // forward (central)
                xhelp[i] = x[i] + dxmin[i];
                dx = xhelp[i]-x[i];
              } else {                        // backwards
                xhelp[i] = x[i] - dxmin[i];
                dx = x[i]-xhelp[i];
              }
              const tdouble tr = dx/x[i];
              if (tr>fd_max) fd_max = tr;
            }
          #if FLX_DEBUG
            if (dx<=ZERO) throw FlxException_Crude("FlxObjFORM_base::do_FORM_6");
          #endif
        if ( RndBox->check_xVec(xhelp) ) {        // check if xhelp is in domain
          RndBox->set_smp_x(xhelp);
          // calc LSF
            zo = LSF->calc()/z_0; ++LSFcalls;
          if (fd_method!=2) {
            if (fd_method==1) {
              dzdy[i] = (zo-z_old)/dx;
            } else {
              dzdy[i] = (z_old-zo)/dx;
            }
          } else {                // central difference
            xhelp *= -ONE;
            xhelp.add(x,2);
            if ( RndBox->check_xVec(xhelp) ) {
              RndBox->set_smp_x(xhelp);
              // calc LSF
                zu = LSF->calc()/z_0; ++LSFcalls;
              // calc dz/dx(i)
                dzdy[i] = (zo-zu)/(2*dx);
            } else {
              dzdy[i] = (zo-z_old)/dx;
            }
          }
        } else {
          if (fd_method!=3) {                // forward (central) -> use backwards
            xhelp[i] = x[i] - dx;
          } else {                        // backwards -> use forwards
            xhelp[i] = x[i] + dx;
          }
          #if FLX_DEBUG
            if ( RndBox->check_xVec(xhelp) == false ) {
              throw FlxException_Crude("FlxObjFORM_base::do_FORM_7");
            }
          #endif
          RndBox->set_smp_x(xhelp);
          zu = LSF->calc()/z_0; ++LSFcalls;
          if (fd_method!=3) {
            dzdy[i] = (z_old-zu)/dx;
          } else {
            dzdy[i] = (zu-z_old)/dx;
          }
        }
      }
      // calculate dz/dy
        dxdy.TransMultVec(dzdy);
    } else {
      throw FlxException_NotImplemented("FlxObjFORM_base::do_FORM_8");
//         // calculate dz/dy
//         for (tuint i = 0; i < DIM; ++i) {
//           xhelp = y;
//           if (std::fabs(y[i]) < fd ) xhelp[i] = y[i] + fd;   // make sure its different from zero
//           else xhelp[i] = (1.0+fd)*y[i];
//           dx = xhelp[i]-y[i];
//           if ( RndBox->check_xVec(xhelp) ) {
//             data->RndBox.set_y_Vec(xhelp);
//             // calc LSF
//               zo = LSF->calc();
//             xhelp = y+y-xhelp;
//             if ( RndBox->check_xVec(xhelp) ) {
//               data->RndBox.set_y_Vec(xhelp);
//               // calc LSF
//                 zu = LSF->calc();
//               // calc dz/dx(i)
//                 dzdy[i] = (zo-zu)/(2.0*dx);
//             } else {
//               dzdy[i] = (zo-z_old)/dx;
//             }
//           } else {
//             xhelp[i] = y[i] - dx;
//             data->RndBox.set_y_Vec(xhelp);
//             zu = LSF->calc();
//             dzdy[i] = (z_old-zu)/dx;
//           }
//         }
    }
    if (only_partial) return dzdy;
    
    if (dzdy.get_NormMax() <= GlobalVar.TOL() ) {
      data->RndCreator.gen_smp(y);
      if (verboseLog) {
        GlobalVar.alert.alert("FlxObjFORM_base::do_FORM_9", "Warning: zero_Gradient -> Initialize random seed." );
      }
    } else {
      dx=((dzdy*y)-z_old)/dzdy.get_Norm2_NOroot();  // this is just a factor
      if (opt_method==1) {
        y = dzdy;
        y *= dx;
      } else {
        d = dzdy;
        d *= dx;
        d -= y;
        // determine c^k
          tdouble c = y.get_Norm2()/dzdy.get_Norm2();
          if (fabs(z_old)>=0.001) {
            xhelp = y; xhelp+=d;
            const tdouble c2 = xhelp.get_Norm2_NOroot()/(2*fabs(z_old));
            if (c2>c) c = c2;
          }
          c *= 2;
        lambda = iHLRF_lambda_start->cast2positive();        // start value for the step-size
        const tdouble epsilon = iHLRF_epsilon->cast2positive();                // parameter of Armijo-rule
        const tdouble reduce = iHLRF_reduce->cast2positive();
        const tdouble m_0 = reduce*y.get_Norm2_NOroot()+c*fabs(z_old);        // current value of merit function
          if (epsilon>=ONE) {
            std::ostringstream ssV;
            ssV << "'iHLRF_epsilon' has to be a value within the interval ]0;1[, and not (" << GlobalVar.Double2String(epsilon) << ")";
            throw FlxException_NeglectInInteractive("FlxObjFORM_base::do_FORM_10", ssV.str() );
          }
        tdouble m_n, m_ne, z_next;                // next value (estimate) of merit-function
        const tuint maxIter_LS = 20;
        tuint i;
        std::stringstream ssV2;
        ssV2 << "epsilon = " << GlobalVar.Double2String(epsilon) << std::endl;
        ssV2 << "m0      = " << GlobalVar.Double2String(m_0) << std::endl;
        ssV2 << "z       = " << GlobalVar.Double2String(z_old) << std::endl;
        for (i=0;i<maxIter_LS;++i) {                // maximum number of iterations
          ssV2 << "iteration: " << i+1 << std::endl;
          ssV2 << "  lambda  = " << GlobalVar.Double2String(lambda) << std::endl;
          // compute estimator
            xhelp = y; 
            xhelp.add(dzdy,c*sign(z_old));
            m_ne = m_0 + epsilon*lambda*(xhelp*d);
          // compute LSF-value
            y_new = y;
            y_new.add(d,lambda);
            RndBox->set_smp(y_new);
            z_next = LSF->calc()/z_0; ++LSFcalls;
            ssV2 << "  z     = " << GlobalVar.Double2String(z_next) << std::endl;
            m_n = (y_new.get_Norm2_NOroot())/2+c*fabs(z_next);
            ssV2 << "m_ne    = " << GlobalVar.Double2String(m_ne) << std::endl;
          if (m_ne > m_n) break;        // get me out of here
          lambda *= reduce;
        }
        if (i==maxIter_LS && m_ne <= m_n) {        // no-convergence
          std::ostringstream ssV;
          ssV << "Maximum number of line-search iterations reached (" << maxIter_LS << ").";
          throw FlxException_NeglectInInteractive("FlxObjFORM_base::do_FORM_11", ssV.str(), ssV2.str() );
        }
        z_old = z_next;
        y = y_new;
      }
    }
    
    beta_new = y.get_Norm2();
    
    // check if we are within reasonable intervals
      bool ok = true;                // false, if y has been reassigned
      dx = rv_Phi(-beta_new);
      if (dx == ZERO) {
        dx = fabs(rv_InvPhi(1e-6));        // this is our new beta
        if (verboseLog) {
          GlobalVar.alert.alert("FlxObjFORM_base::do_FORM_12", "Warning: Interval correction (beta=" + GlobalVar.Double2String(beta_new) + ")" );
        }
        y*=dx/beta_new;
        beta_new = y.get_Norm2();  // equal to: beta_new = dx;
        ok = false;
      }
    
    if (opt_method!=2 || !ok) {
      RndBox->set_smp(y);
    }
    RndBox->get_x_Vec(x.get_tmp_vptr());
    // calc LSF
      if (opt_method!=2 || !ok) {
        z_old = LSF->calc()/z_0;        ++LSFcalls;        // actually this is z_new
      }
      
    // check for maximum number of iterations 
      loopc++;
      if (loopc > maxIter) {
        std::ostringstream ssV;
        ssV << "Maximum number of iterations reached (" << maxIter << ").";
        throw FlxException_NeglectInInteractive("FlxObjFORM_base::do_FORM_14", ssV.str() );
      }
    // compute err_orth
//         err_orth = Norm2(dzdy/Norm2(dzdy)+y/beta_new);
    
    GlobalVar.slogcout(4) << format("  Iter %2i: lsf_err=%9.2e  beta=%6.3f") % loopc % z_old % beta_new;
    GlobalVar.slogcout(4)  << format("  beta_err=%9.2e") %((beta_new-beta_old)/beta_new);
    if (opt_method==2)  {
      GlobalVar.slogcout(4)  << format("  lambda=%8.2e") %lambda;
    }
    GlobalVar.slogcout(4) << std::endl; 
    if (verboseLog) {
      GlobalVar.slog(4) << "    y="; flxVec_simple_plot(GlobalVar.slog(4),y,true,-1,0,true); GlobalVar.slog(4) << std::endl;
      GlobalVar.slog(4) << "    x="; flxVec_simple_plot(GlobalVar.slog(4),x,true,-1,0,true); GlobalVar.slog(4) << std::endl;
    }
  } while ( std::fabs((beta_new-beta_old)/beta_new) > eps1->calc() || fabs(z_old) > eps2->calc() );
      
  if ( dzdy*y > ZERO ) {
    beta_new *= -ONE;
  }
  return dzdy;
}


void FlxObjFORM::task()
{ 
  update_Start();
  // calculate the start vector
    flxVec x(DIM);
    flxVec y(DIM);
    eval_xStart(x);
  // Write log-message
      if (only_partial) {
        GlobalVar.slogcout(4) << "partial_derivative: " << std::endl;
      } else {
        GlobalVar.slogcout(4) << "form: performing FORM analysis. " << std::endl;
      }
  // Perform the FORM interation
    tdouble beta; tuint LSFcalls=0;
    flxVec dzdy(do_FORM(x,y,beta,LSFcalls,only_partial));
  
  // Deal with the results
    if (only_partial) {
      data->ConstMtxBox.insert(cnamey,new FlxSMtx(dzdy));
      GlobalVar.slogcout(3) << "  determined partial derivative:" << std::endl;
      GlobalVar.slog(3) << "    dzdy="; flxVec_simple_plot(GlobalVar.slog(3),dzdy,true,-1,0,true); GlobalVar.slog(3) << std::endl;
      return;
    }
    data->ConstMtxBox.insert(cnamey,new FlxSMtx(y));
    data->ConstMtxBox.insert(cnamex,new FlxSMtx(x));
    const std::string betaDP_constName = betaDP->eval_word(true,true);
    if (betaDP_constName!="") {
      data->ConstantBox.insert(betaDP_constName,beta);
    }
  
  GlobalVar.slogcout(3) << "form: determined design point: " << "betaDP=" << GlobalVar.Double2String(beta) << std::endl;
  if (verboseLog) {
    GlobalVar.slog(3) << "    y="; flxVec_simple_plot(GlobalVar.slog(3),y,true,-1,0,true); GlobalVar.slog(3) << std::endl;
    GlobalVar.slog(3) << "    x="; flxVec_simple_plot(GlobalVar.slog(3),x,true,-1,0,true); GlobalVar.slog(3) << std::endl;
  }
  GlobalVar.slogcout(3) << "  Estimated probability of failure:        " << GlobalVar.Double2String(rv_Phi(-beta)) << std::endl;
  GlobalVar.slogcout(3) << "  Probability of failure 'for sure' within [0; " << GlobalVar.Double2String(ONE-rv_cdf_ChiSquare(DIM,pow2(beta))) << "]" << std::endl;
  GlobalVar.slogcout(3) << "  Total number of LSF-calls:               " << LSFcalls << std::endl;

  sensitivities(y,*RndBox,GlobalVar.slog(3));
}

void FlxObjFORM::sensitivities(const flxVec& y, RBRV_constructor& RndBox, std::ostream& slog, flxVec* svp)
{
  // the original alphas
    if (y.get_N()!=RndBox.get_NRV()) {
      throw FlxException("FlxObjFORM::sensitivities", "Specified vector has wrong dimension.");
    }
    const tdouble beta = y.get_Norm2();
    flxVec sensi = y; sensi/=beta;
  // transformation to correlated space
    flxVec w = sensi;
    RndBox.transform_y2w(sensi.get_tmp_vptr_const(),w.get_tmp_vptr());
    w.normalize();
  // the gamma-approach
    FlxMtxLTri Jinv(y.get_N());
    RndBox.calc_Jinv(Jinv);
  // obtain covariance /hat{x} = hS
      FlxMtx hM(Jinv);
      FlxMtxSym hS(hM.nrows());
      hM.TransposeMmultM(hS);
    // get Jacobian
      Jinv.Invert();
    // multiply with sensitivities
      Jinv.TransMultVec(sensi);
    const tuint N = y.get_N();
    for (tuint i=0;i<N;++i) {
      sensi[i] *= sqrt(hS(i,i));
    }
    sensi.normalize();
    if (svp) {
      *svp = sensi;
    }
 
  slog << "  Sensitivities: \t  gamma\t\t gamma^2\t  %" << std::endl;
  for (tuint i=0;i<N;++i) {
    slog << format("  %-15s\t%9.2e\t%9.2e\t%3.1f") % RndBox.get_rv_name(i) % sensi[i] % pow2(sensi[i]) % (pow2(sensi[i])*100) << "%"  << std::endl;
  }
}

FlxObjFORM::~FlxObjFORM()
{
  delete betaDP;
}

void FlxObjFORM_pdf::task()
{
  FunNumber* fn = new FunNumber(ZERO);
  FunSub* fa = new FunSub(fn,rvfun);
  FlxFunction* lsf = new FlxFunction(fa);
  try {
  LSF = lsf;
  update_Start();
  // calculate the start vector
    flxVec x(DIM);
    flxVec y(DIM);
    eval_xStart(x);
  // evaluate the bounds
    tdouble lbound = lboundF->calc();
    tdouble ubound = uboundF->calc();
    if (ubound<=lbound) {
      std::ostringstream ssV;
      ssV << "The lower bound has to be smaller than the upper bound.";
      throw FlxException("FlxObjFORM_pdf::task_2", ssV.str() );
    }
    tuint ninterval = nintervalF->cast2tuint(false);
    tdouble du = (ubound-lbound)/ninterval;
  // Perform the FORM interation
    tdouble beta;
    tdouble u = lbound;
    tuint LSFcalls = 0;
    while (u<=ubound) {
      fn->set_thenumber(u); 
      flxVec dzdy = do_FORM(x,y,beta,LSFcalls);
      sout() << GlobalVar.Double2String(u,true) << "\t";
      sout() << GlobalVar.Double2String( rv_phi(beta)/dzdy.get_Norm2() ) << "\t";         // PDF
      sout() << GlobalVar.Double2String( rv_Phi(beta) );        //CDF
      sout() << std::endl;      
      u+=du;
    }
  } catch (FlxException &e) {
    FLXMSG("FlxObjFORM_pdf::task_3",1);
    delete fn;
    fa->set_childs_zero();
    delete lsf;
    LSF = NULL;
    throw;
  }
  delete fn;
  fa->set_childs_zero();
  delete lsf;
  LSF = NULL;
}

FlxObjFORM_pdf::FlxObjFORM_pdf(bool dolog, FlxFunction* rvfun, FlxFunction* lboundF, FlxFunction* uboundF, FlxFunction* nintervalF, FlxMtxConstFun* xstart, FlxFunction* fdstepV, FlxFunction* epsdyfF, FlxFunction* eps1, FlxFunction* eps2, FlxFunction* iHLRF_lambda_start, FlxFunction* iHLRF_epsilon, FlxFunction* iHLRF_reduce, unsigned int maxIter, bool verbose, std::string ostreamV, bool dxdyAnalytical, FlxMtxConstFun* dx_min, int fd_method, int opt_method, FlxString* rbrvsets)
: FlxObjFORM_base(dolog,NULL,fdstepV,epsdyfF,eps1,eps2,iHLRF_lambda_start,iHLRF_epsilon,iHLRF_reduce,maxIter,verbose,ostreamV,dxdyAnalytical,xstart,dx_min,fd_method,opt_method,rbrvsets), rvfun(rvfun), lboundF(lboundF), uboundF(uboundF), nintervalF(nintervalF)
{
  
}

FlxObjFORM_pdf::~FlxObjFORM_pdf()
{
  delete rvfun;
  delete lboundF;
  delete uboundF;
  delete nintervalF;
}

FlxObjReadKDE::FlxObjReadKDE()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(false,"kde::do_cdf"));
  ParaBox.insert("do_cdf", "kde::do_cdf" );
  AllDefParaBox->insert(new FlxOptionalParaFun(-10,"kde::lbound"));
  ParaBox.insert("lbound", "kde::lbound" );
  AllDefParaBox->insert(new FlxOptionalParaFun(10,"kde::ubound"));
  ParaBox.insert("ubound", "kde::ubound" );
  AllDefParaBox->insert(new FlxOptionalParaFun(100,"kde::ninterval"));
  ParaBox.insert("ninterval", "kde::ninterval" );
}

FlxObjBase* FlxObjReadKDE::read()
{
  FlxFunction* funR = NULL; FlxString* rbrvsets = NULL; FlxFunction* N=NULL; FlxRndKernel_base* kernel=NULL; FlxFunction* h=NULL;
  try {
    reader->getChar('(',false);
    funR = new FlxFunction(funReader,false);
    reader->getChar(',',false);
    rbrvsets = new FlxString(false,false);
    reader->getChar(',',false);
    N = new FlxFunction(funReader,false);
    reader->getChar(',',false);
    kernel = FlxRndKernel_base::read(false);
    reader->getChar(',',false);
    h = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjKDE(get_doLog(),funR,rbrvsets, N, kernel, h, get_optPara_FlxFunction("lbound"), get_optPara_FlxFunction("ubound"), get_optPara_FlxFunction("ninterval"), get_optPara_bool("do_cdf"), get_stream());
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadKDE::read",1);
    if (funR) delete funR;
    if (rbrvsets) delete rbrvsets;
    if (N) delete N;
    if (kernel) delete kernel;
    if (h) delete h;
    throw;
  }
}

FlxObjKDE::FlxObjKDE( bool dolog, FlxFunction* funR, FlxString* rbrvsets, FlxFunction* N, FlxRndKernel_base* kernel, FlxFunction* h, FlxFunction* lbound, FlxFunction* ubound, FlxFunction* Ninterval, bool do_cdf, std::string ostreamV) 
: FlxObjOutputBase(dolog,ostreamV), funR(funR), rbrvsets(rbrvsets), N(N), kernel(kernel), h(h), lbound(lbound), ubound(ubound), Ninterval(Ninterval), do_cdf(do_cdf)
{
  
}

FlxObjKDE::~FlxObjKDE()
{
  delete funR;
  delete rbrvsets;
  delete N;
  delete kernel;
  delete h;
  delete lbound;
  delete ubound;
  delete Ninterval;
}

void FlxObjKDE::task()
{
  const std::string setstr = rbrvsets->eval(true);
  const std::vector<std::string> set_str_vec = parse_strseq_as_vec(setstr);
  RBRV_constructor RndBox(set_str_vec,data->rbrv_box);
  
  tdouble dlb = lbound->calc();          // lower bound
  tdouble dub = ubound->calc();                // upper bound
  if (dlb >= dub) {
    std::ostringstream ssV_2;
    ssV_2 << "Lower bound (" << dlb << ") must not be larger than upper bound (" << dub << ").";
    throw FlxException_NeglectInInteractive("FlxObjKDE::task_1", ssV_2.str() ); 
  }
  tulong LNi = Ninterval->cast2tuint(false);                 // number of intervals
  #if FLX_KAHAN_KDE
    pdouble dss(dub);
    dss-=dlb;
    dss/=tdouble(LNi);
  #else
    tdouble dss = (dub-dlb)/LNi;                                 // stepsize
  #endif
  kernel->set_h(h->cast2positive(false));
  tulong LN = N->cast2tulong();                                // number of samples
  
  // generate random realizations
  tdouble *arr = new tdouble[LN];                        // array
  for (tulong i = 0; i < LN; ++i) {
    RndBox.gen_smp();
    arr[i] = funR->calc();
  }
  
  #if FLX_KAHAN_KDE
    pdouble x = dlb;
    pdouble p;
  #else
    tdouble x = dlb;
    tdouble p;
  #endif
  if (!do_cdf) {
    for (tulong i = 0; i <= LNi; ++i) {
      p = ZERO;
      for (tulong j = 0; j < LN; ++j) {
        #if FLX_KAHAN_KDE
          p+=kernel->pdf(x.cast2double(), arr[j]);
        #else
          p+=kernel->pdf(x, arr[j]);
        #endif
      }
      p/=tdouble(LN);
      #if FLX_KAHAN_KDE
        sout() << GlobalVar.Double2String(x.cast2double()) << " " << GlobalVar.Double2String(p.cast2double()) << std::endl;
      #else
        sout() << GlobalVar.Double2String(x) << " " << GlobalVar.Double2String(p) << std::endl;
      #endif
      x += dss;
    }
  } else  {
    for (tulong i = 0; i <= LNi; ++i) {
      p = ZERO;
      for (tulong j = 0; j < LN; ++j) {
        #if FLX_KAHAN_KDE
          p+=kernel->cdf(x.cast2double(), arr[j]);
        #else
          p+=kernel->cdf(x, arr[j]);
        #endif
      }
      p/=tdouble(LN);
      #if FLX_KAHAN_KDE
        sout() << GlobalVar.Double2String(x.cast2double()) << " " << GlobalVar.Double2String(p.cast2double()) << std::endl;
      #else
        sout() << GlobalVar.Double2String(x) << " " << GlobalVar.Double2String(p) << std::endl;
      #endif
      x += dss;
    }
  }
  
  delete[] arr;
  
  GlobalVar.slog(4) << "kde: performed kernel density estimation (" << (do_cdf?"cdf":"pdf") << ") of expression '" << funR->write() << "' (" << setstr << ")." << std::endl;
  GlobalVar.slog(4) << "  Kernel: "; kernel->print_info(GlobalVar.slog(4)); GlobalVar.slog(4) << " (Number of samples: " << GlobalVar.Double2String(LN) << ")" << std::endl;
  GlobalVar.slog(4) << "  Plot interval: [" << GlobalVar.Double2String(dlb) << "; " << GlobalVar.Double2String(dub) << "]; stepsize=" 
    #if FLX_KAHAN_KDE
      << dss.cast2double() << std::endl;
    #else
      << dss << std::endl;
    #endif
}


FlxObjReadMCSsensi::FlxObjReadMCSsensi(): FlxObjReadOutputBase()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(1e2,"mcssensi::nb"));
  ParaBox.insert("nb", "mcssensi::nb" );
}

FlxObjBase* FlxObjReadMCSsensi::read()
{
  FlxMtxConstFun* mcn = new FlxMtxConstFun(false);
  FlxString* istrm = NULL;
  FlxFunction* rvN = NULL;
  try {
    reader->getChar('=');
    istrm = new FlxString(false,false);
    reader->getChar('(');
    rvN = new FlxFunction(funReader,false);
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjMCSsensi(get_doLog(), get_stream(),mcn,istrm,rvN,get_optPara_FlxFunction("nb"));
  } catch (FlxException& e) {
    delete mcn;
    if (istrm) delete istrm;
    if (rvN) delete rvN;
    throw;
  }
}

FlxObjMCSsensi::FlxObjMCSsensi(bool dolog, std::string ostreamV, FlxMtxConstFun* mcn, FlxString* istrm, FlxFunction* rvN, FlxFunction* nb)
: FlxObjOutputBase(dolog,ostreamV), mcn(mcn), istrm(istrm), rvN(rvN), nb(nb)
{

}

FlxObjMCSsensi::~FlxObjMCSsensi()
{
  delete mcn;
  delete istrm;
  delete rvN;
  delete nb;
}

void FlxObjMCSsensi::task()
{
  const tuint M = rvN->cast2tuint(false);                // number of random variables
  const tuint Lb = nb->cast2tuint(false);                // length of blocks
  const std::string isn = istrm->eval_word(true);
  FlxIstream& is = data->IstreamBox.get(isn);
  FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&is);
  if (isv==NULL) {
    throw FlxException_NeglectInInteractive("FlxObjMCSsensi::task_01", "Input stream '" + isn + "' is not a vector-input stream.");
  }
  // use vector 'tv' to point to entire data-set
    const tdouble* tv = isv->get_tmpPtr();
    const size_t N = isv->get_total_size();                // total number of entries
    if ( N % (M+1) != 0) {
      throw FlxException_NeglectInInteractive("FlxObjMCSsensi::task_02", "There appears to be an issue with the size of the input stream.");
    }
    const size_t Ne = N / (M+1);                        // total number of entries
    const tuint Nb = Ne/Lb;                                // number of blocks
    tuint Nb_rem = Ne % Lb;                        // division remainder for distribution in boxes
  // add limit-state values to 'vec_lsf'
    flxVec vec_lsf(Ne);
    const tdouble* tv1 = tv;
    for (size_t i=0;i<Ne;++i) {
      vec_lsf[i] = *tv1;
      tv1 += (M+1);
    }
    const tdouble lsf_mu = vec_lsf.get_Mean();
    const tdouble lsf_var = vec_lsf.get_Var(lsf_mu);
    flxVec sensi_vec(M);
  std::vector<kvEl> sbox;
  kvEl el;
  for (tuint i=0;i<M;++i) {
    // sort with respect to random variable 'i'
      sbox.clear();
      tv1 = tv;
      for (size_t j=0;j<Ne;++j) {
          const tdouble lsf_val = *tv1;
          #if FLX_DEBUG
            if (vec_lsf[j] != lsf_val) throw FlxException_Crude("FlxObjMCSsensi::task_03");
          #endif
          tv1 += i+1;
          const tdouble rv_val = *tv1;
          tv1 += size_t(M-i);
          el.key = rv_val;
          el.val = lsf_val;
          sbox.push_back(el);
      }
      std::sort(sbox.begin(), sbox.end());
    // compute block mean and variance 
      size_t j = 0;
      flxVec bmv(Nb);
      flxVec bvv(Nb);
      for (size_t k=0;k<Nb;++k) {
        size_t Lb_t = Lb;
        if (Nb_rem>0) {
          ++Lb_t;
          --Nb_rem;
        }
        flxVec tlv(Lb_t);
        for (size_t m=0;m<Lb_t;++m) {
          #if FLX_DEBUG
            if (j>=Ne) {
              throw FlxException_Crude("FlxObjMCSsensi::task_04");
            }
          #endif
          tlv[m] = sbox[j++].val;
        }
        bmv[k] = tlv.get_Mean();
        bvv[k] = tlv.get_Var(bmv[k]);
      }
      const tdouble bmv_mean = bmv.get_Mean();
      const tdouble bmv_var = bmv.get_Var(bmv_mean);
      const tdouble bvv_mean = bvv.get_Mean();
      const tdouble bmv_var_corrected = bmv_var * (lsf_var/(bmv_var+bvv_mean));
      sout() << GlobalVar.Double2String(i) << " " << GlobalVar.Double2String(lsf_var/(bmv_var+bvv_mean)) << std::endl;
      sensi_vec[i] = bmv_var_corrected/lsf_var;
      
  }
  
  tuint M_ = M;
  tdouble* stvp = data->ConstMtxBox.get_Vec(mcn->eval(),M_);
  flxVec stv(stvp,M);
  stv = sensi_vec;
  
  sout() << "first-order sensitivities (sum=" << GlobalVar.Double2String(sensi_vec.get_sum()) << "):" << std::endl;
  sout() << sensi_vec << std::endl;
}


FlxObjReadStatSmp::FlxObjReadStatSmp()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(1e6,"statsmp::np"));
  ParaBox.insert("np", "statsmp::np" );
  
  AllDefParaBox->insert(new FlxOptionalParaFlxString("","statsmp::addname",false));
  ParaBox.insert("addname", "statsmp::addname" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"statsmp::optionp"));
  ParaBox.insert("optionp", "statsmp::optionp" );
  
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"statsmp::sigfig"));
  ParaBox.insert("sigfig", "statsmp::sigfig" );
}

FlxObjBase* FlxObjReadStatSmp::read()
{
  FlxString* flxStr = NULL;
  FlxString* aStr = NULL;
  try {
    reader->getChar('(',false);
    flxStr = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    const int optionp = get_optPara_int_from_FlxFunction("optionp");
    aStr = get_optPara_FlxString("addname");
    // define constants    
      if (flxStr->is_static()) {
        const std::string as = aStr->eval_word(true,true);
        if (as.empty()==false) {
          data->ConstantBox.declareC(as + "_n");
          data->ConstantBox.declareC(as + "_mean");
          data->ConstantBox.declareC(as + "_mean_cov");
          data->ConstantBox.declareC(as + "_sd");
          data->ConstantBox.declareC(as + "_coeffofvar");
          data->ConstantBox.declareC(as + "_var");
          data->ConstantBox.declareC(as + "_skewness");
          data->ConstantBox.declareC(as + "_min");
          data->ConstantBox.declareC(as + "_max");
          data->ConstantBox.declareC(as + "_range");
          if (optionp>1) {
            data->ConstantBox.declareC(as + "_mean2");
            data->ConstantBox.declareC(as + "_mean2_cov");
            data->ConstantBox.declareC(as + "_sd2");
            data->ConstantBox.declareC(as + "_coeffofvar2");
            data->ConstantBox.declareC(as + "_var2");
            data->ConstantBox.declareC(as + "_skewness2");
            data->ConstantBox.declareC(as + "_min2");
            data->ConstantBox.declareC(as + "_max2");
            data->ConstantBox.declareC(as + "_range2");
            data->ConstantBox.declareC(as + "_corrcoeff");
            data->ConstantBox.declareC(as + "_covar");
          }
        }
      }
    return new FlxObjStatSmp(get_doLog(), get_stream(), flxStr, aStr ,get_optPara_FlxFunction("np"), optionp, get_optPara_bool("sigfig"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadStatSmp::read",1);
    if (flxStr) delete flxStr;
    if (aStr) delete aStr;
    throw;
  }
}


struct sigfig_mean_rss {
  tdouble mu;
  tdouble sd;
  tdouble ubound;
  tdouble lbound;
  tdouble p;
};

tdouble sigfig_mean_root_search_fun(const tdouble N, void* dp) {
  sigfig_mean_rss *s = (sigfig_mean_rss *)dp;
  const tdouble f = sqrt(N)/s->sd;
  if (N<=1.5) return s->p - ZERO;
  else if (std::isinf(N)) return s->p - ONE;
  return s->p - (rv_cdf_Studentst(N-1,(s->ubound-s->mu)*f) - rv_cdf_Studentst(N-1,(s->lbound-s->mu)*f));
}

void FlxObjStatSmp::sigfig_mean(const tdouble estim_mu, const tdouble estim_sd, const tulong N)
{
  if (fabs(estim_mu)<=GlobalVar.TOL()) return;
  if (estim_sd<=GlobalVar.TOL()) return;
  const tdouble f = sqrt(N)/estim_sd;
  sout() << "  Pr[significant figures of mean]:         " << std::endl;
  int p10 = std::log10(std::fabs(estim_mu));
    if (fabs(estim_mu)<ONE) --p10;
  for (int i=0;i<=5;++i) {
    const tdouble rd = round_flx_fb(estim_mu,i);
    #ifdef __unix__
      const tdouble rd_delta_h = exp10(p10-i)/2;
    #else
      const tdouble rd_delta_h = pow(10,p10-i)/2;
    #endif
    const tdouble Pr_sf = rv_cdf_Studentst(tdouble(N-1),((rd-estim_mu)+rd_delta_h)*f) - rv_cdf_Studentst(tdouble(N-1),((rd-estim_mu)-rd_delta_h)*f);
    sout() << "    " << GlobalVar.Double2String(i,false,-1,5) << GlobalVar.Double2String_sci(rd,i,20) << "      " << GlobalVar.Double2String(Pr_sf) << std::endl;
  } 
  sout() << "  Nsamples needed to maintain specified accuracy of mean:         " << std::endl;
  sout() << "            accuracy                 50%        75%        90%        95%        99% " << std::endl;
  sigfig_mean_rss ts;
  ts.mu = estim_mu;
  ts.sd = estim_sd;
  for (int i=0;i<=5;++i) {
    const tdouble rd = round_flx_fb(estim_mu,i);
    #ifdef __unix__
      const tdouble rd_delta_h = exp10(p10-i)/2;
    #else
      const tdouble rd_delta_h = pow(10,p10-i)/2;
    #endif
    ts.ubound = rd + rd_delta_h;
    ts.lbound = rd - rd_delta_h;
    sout() << "    " << GlobalVar.Double2String(i,false,-1,5) << GlobalVar.Double2String_sci(rd,i,20) << "      ";
    ts.p = 0.5;
    sout() << GlobalVar.Double2String(round_flx(flx_RootSearch_RegulaFalsi(&sigfig_mean_root_search_fun,&ts,tdouble(10),tdouble(1e6))),false,3,10) << " ";
    ts.p = 0.75;
    sout() << GlobalVar.Double2String(round_flx(flx_RootSearch_RegulaFalsi(&sigfig_mean_root_search_fun,&ts,tdouble(10),tdouble(1e6))),false,3,10) << " ";
    ts.p = 0.9;
    sout() << GlobalVar.Double2String(round_flx(flx_RootSearch_RegulaFalsi(&sigfig_mean_root_search_fun,&ts,tdouble(10),tdouble(1e6))),false,3,10) << " ";
    ts.p = 0.95;
    sout() << GlobalVar.Double2String(round_flx(flx_RootSearch_RegulaFalsi(&sigfig_mean_root_search_fun,&ts,tdouble(10),tdouble(1e6))),false,3,10) << " ";
    ts.p = 0.99;
    sout() << GlobalVar.Double2String(round_flx(flx_RootSearch_RegulaFalsi(&sigfig_mean_root_search_fun,&ts,tdouble(10),tdouble(1e6))),false,3,10) << " ";
    sout() << std::endl;
  }
}

void FlxObjStatSmp::task()
{
  bool doubleSample;
  if (optionP==1) doubleSample = false;
  else if (optionP==2) doubleSample = true;
  else {
    throw FlxException_Crude("FlxObjStatSmp::task_1");
  }
  
  const std::string as = (add_name?(add_name->eval_word(true,true)):add_name_str);
  FlxIstream &is = isname?(data->IstreamBox.get(isname->eval_word(true))):(*is_vec);
  tdouble d,d2;
  tulong N=0;
  tuint NpB = tuint(sqrt(tdouble(Np?(Np->cast2tuint(false)):Np_val)));
  qdouble sx1(NpB,true), sx1_2(NpB,true);
  qdouble sx2(NpB,true), sx2_2(NpB,true);
  qdouble sx3(NpB,true), sx3_2(NpB,true);
  qdouble sxy(NpB,true);
  tdouble xmin,xmax,xmin2,xmax2;
  if (is.get_value(d,true)) {
    xmax = d; xmin = d;
    if (doubleSample) {
      if (!is.get_value(d2,true)) {
        std::ostringstream ssV_2;
        ssV_2 << "The two sets have not the same size.";
        throw FlxException_NeglectInInteractive("FlxObjStatSmp::task_2", ssV_2.str() );
      }
      xmax2 = d2; xmin2 = d2;
    }
    do {
      ++N;
      sx1+=d;
      sx2+=pow2(d);
      sx3+=pow_int(d,3);
      if (d<xmin) xmin=d;
      if (d>xmax) xmax=d;
      if (doubleSample) {
        if (N>1) {
          if (!is.get_value(d2,true)) {
            std::ostringstream ssV_2;
            ssV_2 << "The two sets have not the same size.";
            throw FlxException_NeglectInInteractive("FlxObjStatSmp::task_3", ssV_2.str() );
          }
        }
        sx1_2+=d2;
        sx2_2+=pow2(d2);
        sx3_2+=pow_int(d2,3);
        sxy+=d*d2;
        if (d2<xmin2) xmin2=d2;
        if (d2>xmax2) xmax2=d2;
      }
    } while (is.get_value(d,true));
  } else {
    std::ostringstream ssV_2;
    ssV_2 << "No sample points.";
    throw FlxException_NeglectInInteractive("FlxObjStatSmp::task_4", ssV_2.str() );
  }
  // mean
    pdouble t = sx1.cast2pdouble();
    t/=N;
    pdouble mean = t;
    pdouble mean2, t_2;
    if (doubleSample) {
      t_2 = sx1_2.cast2pdouble();
      t_2/=N;
      mean2 = t_2;
    }
  // variance
    t*=mean;
    pdouble t2 = sx2.cast2pdouble();
    t2-=t*N;
    if (t2.cast2double()<ZERO) {
      if (fabs(t2.cast2double())>GlobalVar.TOL()) GlobalVar.alert.alert("FlxObjStatSmp::task_5","a negative variance should not occur");
      t2 = ZERO;
    }
    pdouble variance = t2;
    variance/=(N-ONE);
    pdouble variance2;
    if (doubleSample) {
      t_2*=mean2;
      t2 = sx2_2.cast2pdouble();
      t2-=t_2*N;
      if (t2.cast2double()<ZERO) {
        if (fabs(t2.cast2double())>GlobalVar.TOL()) GlobalVar.alert.alert("FlxObjStatSmp::task_6","a negative variance should not occur");
        t2 = ZERO;
      }
      variance2 = t2;
      variance2/=(N-ONE);
    }
  // skewness
    t2 = t;
    t2*=mean;
    t2*=(2*ONE)*N;
    t2+=sx3.cast2pdouble();
    t=sx2.cast2pdouble();
    t*=mean;
    t*=tdouble(3.0);
    t2-=t;
    t=variance;
    t*=N*sqrt(variance.cast2double());
    tdouble skewness = t2.cast2double()/t.cast2double();
    tdouble skewness2;
    if (doubleSample) {
      t2 = t_2;
      t2*=mean2;
      t2*=(2*ONE)*N;
      t2+=sx3_2.cast2pdouble();
      t_2=sx2_2.cast2pdouble();
      t_2*=mean2;
      t_2*=tdouble(3.0);
      t2-=t_2;
      t_2=variance2;
      t_2*=N*sqrt(variance2.cast2double());
      skewness2 = t2.cast2double()/t_2.cast2double();
    }
    tdouble covariance,corrcoeff;
    if (doubleSample) {
      t2 = sxy.cast2pdouble();
      t = mean;
      t *= mean2;
      t *= N;
      t2 -= t;
      covariance = t2.cast2double() / (N-ONE);
      corrcoeff = covariance / sqrt( variance.cast2double() * variance2.cast2double() );
    }
    // estim. coefficient of variation mean
      const tdouble mean_cov = sqrt(variance.cast2double())/(sqrt(tdouble(N))*fabs(mean.cast2double()));
      const tdouble mean2_cov = sqrt(variance2.cast2double())/(sqrt(tdouble(N))*fabs(mean2.cast2double()));
  if (NOTdolog==false) {
    sout() << "Evaluation of the statistical properties of the given samples:" << std::endl;
    sout() << "  number of samples:               " << GlobalVar.Double2String(N) << std::endl;
    if (doubleSample) {
    sout() << std::endl;
    sout() << "  Sample Set 1" << std::endl;
    sout() << "  ------------" << std::endl;
    }
    sout() << "  sample mean:                     " << GlobalVar.Double2String(mean.cast2double()) << "\t(estim. c.o.v.: " << GlobalVar.Double2String(mean_cov) << ")" << std::endl;
    sout() << "  sample standard deviation:       " << GlobalVar.Double2String(sqrt(variance.cast2double())) << std::endl;
    if (mean.cast2double()!=ZERO) {
    sout() << "  sample coefficient of variation: " << GlobalVar.Double2String(sqrt(variance.cast2double())/fabs(mean.cast2double())) << std::endl;
    }
    sout() << "  sample variance:                 " << GlobalVar.Double2String(variance.cast2double()) << std::endl;
    sout() << "  sample skewness:                 " << GlobalVar.Double2String(skewness) << std::endl;
    sout() << "  smallest sample:                 " << GlobalVar.Double2String(xmin) << std::endl;
    sout() << "  largest sample:                  " << GlobalVar.Double2String(xmax) << std::endl;
    sout() << "  sample range:                    " << GlobalVar.Double2String(xmax-xmin) << std::endl;
    if (sigfig) {
      sigfig_mean(mean.cast2double(),sqrt(variance.cast2double()),N);
    }
    if (doubleSample) {    
    sout() << std::endl;
    sout() << "  Sample Set 2" << std::endl;
    sout() << "  ------------" << std::endl;
    sout() << "  sample mean:                     " << GlobalVar.Double2String(mean2.cast2double()) << "\t(estim. c.o.v.: " << GlobalVar.Double2String(mean2_cov) << ")" << std::endl;
    sout() << "  sample standard deviation:       " << GlobalVar.Double2String(sqrt(variance2.cast2double())) << std::endl;
    if (mean2.cast2double()!=ZERO) {
    sout() << "  sample coefficient of variation: " << GlobalVar.Double2String(sqrt(variance2.cast2double())/fabs(mean2.cast2double())) << std::endl;
    }
    sout() << "  sample variance:                 " << GlobalVar.Double2String(variance2.cast2double()) << std::endl;
    sout() << "  sample skewness:                 " << GlobalVar.Double2String(skewness2) << std::endl;
    sout() << "  smallest sample:                 " << GlobalVar.Double2String(xmin2) << std::endl;
    sout() << "  largest sample:                  " << GlobalVar.Double2String(xmax2) << std::endl;
    sout() << "  sample range:                    " << GlobalVar.Double2String(xmax2-xmin2) << std::endl;
    sout() << std::endl;
    sout() << "  Sample Dependence" << std::endl;
    sout() << "  -----------------" << std::endl;
    sout() << "  sample correlation coefficient:  " << GlobalVar.Double2String(corrcoeff) << std::endl;
    sout() << "  sample covariance:               " << GlobalVar.Double2String(covariance) << std::endl;
    }
  }
  if (!as.empty()) {
    data->ConstantBox.insert(as + "_n",N);
    data->ConstantBox.insert(as + "_mean",mean.cast2double());
    data->ConstantBox.insert(as + "_mean_cov",mean_cov);
    data->ConstantBox.insert(as + "_sd",sqrt(variance.cast2double()));
    if (mean.cast2double()!=ZERO) {
    data->ConstantBox.insert(as + "_coeffofvar",sqrt(variance.cast2double())/fabs(mean.cast2double()));
    }
    data->ConstantBox.insert(as + "_var",variance.cast2double());
    data->ConstantBox.insert(as + "_skewness",skewness);
    data->ConstantBox.insert(as + "_min",xmin);
    data->ConstantBox.insert(as + "_max",xmax);
    data->ConstantBox.insert(as + "_range",xmax-xmin);
    if (doubleSample) {
      data->ConstantBox.insert(as + "_mean2",mean2.cast2double());
    data->ConstantBox.insert(as + "_mean2_cov",mean2_cov);
      data->ConstantBox.insert(as + "_sd2",sqrt(variance2.cast2double()));
      if (mean2.cast2double()!=ZERO) {
      data->ConstantBox.insert(as + "_coeffofvar2",sqrt(variance2.cast2double())/fabs(mean2.cast2double()));
      }
      data->ConstantBox.insert(as + "_var2",variance2.cast2double());
      data->ConstantBox.insert(as + "_skewness2",skewness2);
      data->ConstantBox.insert(as + "_min2",xmin2);
      data->ConstantBox.insert(as + "_max2",xmax2);
      data->ConstantBox.insert(as + "_range2",xmax2-xmin2);
      data->ConstantBox.insert(as + "_corrcoeff",corrcoeff);
      data->ConstantBox.insert(as + "_covar",covariance);
    }
  }
}

FlxObjStatSmp::FlxObjStatSmp(bool dolog, std::string ostreamV, FlxString* isname, FlxString* add_name, FlxFunction* Np, int optionP, const bool sigfig)
: FlxObjOutputBase(dolog,ostreamV), isname(isname), is_vec(NULL), add_name(add_name), add_name_str(""), Np(Np), Np_val(0), optionP(optionP), sigfig(sigfig)
{

}

FlxObjStatSmp::FlxObjStatSmp(bool dolog, std::string ostreamV, FlxIstream_vector* is_vec, const std::string& add_name_str, const tuint Np_val, int optionP, const bool sigfig)
: FlxObjOutputBase(dolog,ostreamV), isname(NULL), is_vec(is_vec), add_name(NULL), add_name_str(add_name_str), Np(NULL), Np_val(Np_val), optionP(optionP), sigfig(sigfig) 
{
  
}

FlxObjStatSmp::~FlxObjStatSmp()
{
  if (isname) delete isname;
  if (add_name) delete add_name;
  if (Np) delete Np;
}

FlxObjReadSortSmp::FlxObjReadSortSmp()
{
  AllDefParaBox->insert(new FlxOptionalParaFun(1e6,"sortsmp::np"));
  ParaBox.insert("np", "sortsmp::np" );
}

FlxObjBase* FlxObjReadSortSmp::read()
{
  FlxString* flxStr = NULL;
  try {
    reader->getChar('(',false);
    flxStr = new FlxString(false,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjSortSmp(get_doLog(), get_stream(), flxStr, get_optPara_FlxFunction("np"));
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSortSmp::read",1);
    if (flxStr) delete flxStr;
    throw;
  }
}

FlxObjSortSmp::FlxObjSortSmp(bool dolog, std::string ostreamV, FlxString* isname, FlxFunction* Np)
: FlxObjOutputBase(dolog,ostreamV), isname(isname), Np(Np)
{
  
}

FlxObjSortSmp::~FlxObjSortSmp()
{
  delete isname;
  delete Np;
}

void FlxObjSortSmp::task()
{
  const std::string isn = isname->eval_word(true);
  FlxIstream &is = data->IstreamBox.get(isn);
  FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&is);
  if (isv!=NULL) {
    isv->sortStream();
    GlobalVar.slog(3) << "sortsmp: sorted vector-input stream '" << isn << "' with " << GlobalVar.Double2String(isv->get_total_size()) << " entries." << std::endl;
    isv->reset_stream();
  } else {
    const tulong BS = Np->cast2tulong(false);
    
    std::vector<tdouble> sortVec;
    sortVec.reserve(BS);
    tdouble d;
    while (is.get_value(d,true)) {
      sortVec.push_back(d);
    };
    std::sort(sortVec.begin(),sortVec.end());
    std::vector<tdouble>::iterator it;
    for (it=sortVec.begin(); it!=sortVec.end(); ++it) {
      sout() << GlobalVar.Double2String(*it) << std::endl;
    }
  }
}


FlxObjReadSmpPlot::FlxObjReadSmpPlot()
{
  AllDefParaBox->insert(new FlxOptionalParaBool(true,"smpplot::autobound"));
  ParaBox.insert("autobound", "smpplot::autobound" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"smpplot::xmin"));
  ParaBox.insert("xmin", "smpplot::xmin" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"smpplot::xmax"));
  ParaBox.insert("xmax", "smpplot::xmax" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"smpplot::binestimator"));
  ParaBox.insert("binestimator", "smpplot::binestimator" );
  
  AllDefParaBox->insert(new FlxOptionalParaFun(ZERO,"smpplot::nbins"));
  ParaBox.insert("nbins", "smpplot::nbins" );
}

FlxObjBase* FlxObjReadSmpPlot::read()
{
  FlxString* flxStr = NULL;
  FlxFunction* typeFun = NULL;
  try {
    reader->getChar('(',false);
    flxStr = new FlxString(false,false);
    reader->getChar(',',false);
    reader->getWord("type",false);
    reader->getChar('=',false);
    typeFun = new FlxFunction(funReader,false);
    reader->getChar(')',false);
    read_optionalPara(false);
    const int binEstimator = get_optPara_int_from_FlxFunction("binestimator");
    return new FlxObjSmpPlot(get_doLog(),get_stream(),flxStr,typeFun,get_optPara_bool("autobound"),get_optPara_FlxFunction("xmin"),get_optPara_FlxFunction("xmax"),binEstimator,get_optPara_FlxFunction("nbins"),get_prec(),get_fixW());
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadSmpPlot::read",1);
    if (flxStr) delete flxStr;
    if (typeFun) delete typeFun;
    throw;
  }
}

FlxObjSmpPlot::FlxObjSmpPlot(const bool dolog, const std::string ostreamV, FlxString* isname, FlxFunction* typeFun, const bool autoBound, FlxFunction* xminU, FlxFunction* xmaxU, const int binEstimator, FlxFunction* NbinsU, const int prec, const int fixW)
: FlxObjOutputBase(dolog, ostreamV,false,false,prec,fixW), isname(isname), typeFun(typeFun), autoBound(autoBound), xminU(xminU), xmaxU(xmaxU), binEstimator(binEstimator), NbinsU(NbinsU)
{

}

FlxObjSmpPlot::~FlxObjSmpPlot()
{
  delete isname;
  delete typeFun;
  delete xminU;
  delete xmaxU;
  delete NbinsU;
}

void FlxObjSmpPlot::task()
{
  const int type = typeFun->cast2int();
  const std::string isname_str = isname->eval_word(true);
  FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(isname_str)));
  if (isv==NULL) {
    std::ostringstream ssV;
    ssV << "The input stream'" << isname_str << "' is not a vector-input stream.";
    throw FlxException_NeglectInInteractive("FlxObjSmpPlot::task_1", ssV.str() );
  }
  isv->reset_stream();
  const tulong N = isv->get_total_size();
  GlobalVar.slog(4) << "smpplot: vector-input stream with " << GlobalVar.Double2String(isv->get_total_size()) << " entries." << std::endl;
  
  if (type==1) {        // cumulative frequency diagram
    tdouble pdxl = -log(tdouble(N));
    tulong count = 0;
    tdouble px = ZERO;
    tdouble qpx = ONE;
    tdouble d;
    isv->get_value(d,true);
    tdouble d_prev = d;
    bool b = true;
    while (b) {
      d_prev = d;
      sout() << GlobalVar.Double2String(d,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(px,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(qpx,checkTOL,prec,fixW) << std::endl;
      while (d==d_prev) {
        ++count;
        px = exp(pdxl + log(tdouble(count)));
        qpx = exp(pdxl + log(tdouble(N-count)));
        if ((b=isv->get_value(d,true))==false) break;
      }
      if (!b) qpx = ZERO;
      sout() << GlobalVar.Double2String(d_prev,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(px,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(qpx,checkTOL,prec,fixW) << std::endl;
    }
  } else if (type==2 || type==3 || type==4 || type==5 ) {        // histogram
    // get the bounds
      tdouble min, max;
      if (autoBound) {
        min=isv->get_first_value();
        max=isv->get_last_value();
        if (min>=max) {
        std::ostringstream ssV;
        ssV << "Set vector-input stream appears not to be sorted.";
        throw FlxException_NeglectInInteractive("FlxObjSmpPlot::task_3", ssV.str() );
      }
      } else {
        min=xminU->calc();
        max=xmaxU->calc();
        if (min>=max) {
        std::ostringstream ssV;
        ssV << "xmin(" << min << ") has to be smaller than xmax(" << max << ").";
        throw FlxException_NeglectInInteractive("FlxObjSmpPlot::task_4", ssV.str() );
      }
      }
    // get the number of bins
      tuint Nbins = NbinsU->cast2tuintW0(false);
      if (Nbins==0) {
        if (binEstimator==1) {
          Nbins = tuint(round_flx( sqrt((tdouble)N) )) ;
        } else if (binEstimator==2) {
          Nbins = tuint(round_flx( ONE + 3.3 * log10(tdouble(N)) )) ;
        } else {
          std::ostringstream ssV;
          ssV << "'binEstimator' has to be eigher 1 or 2 and not '" << binEstimator << "'.";
          throw FlxException_NeglectInInteractive("FlxObjSmpPlot::task_5", ssV.str() );
        }
        if (Nbins<5) Nbins = 5;
        else if (Nbins > 25) Nbins = 25;
      }
    const tdouble dx = (max-min)/Nbins;
    tdouble xm = min + dx/2;
    tdouble xe = min + dx;
    tulong xn = 0;
    
    tdouble d;
    bool stream_empty = false;
    if (!isv->get_value(d,true)) stream_empty = true;
    if (type == 5 ) {
      sout() << GlobalVar.Double2String(xm-dx,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(ZERO,checkTOL,prec,fixW) << std::endl;
    }
    for (tuint i=1;i<=Nbins;++i) {
      do {
        if (!stream_empty) {
          if (d<xe) {
            ++xn;
          } else {
            break;
          }
        } else break;
        if (!isv->get_value(d,true)) stream_empty = true;
      } while (true);
      if (i==Nbins && !stream_empty) {
        while (d<=xe) {
          ++xn;
          if (!isv->get_value(d,true)) {
            stream_empty = true;
            break;
          }
        }
      }
      sout() << GlobalVar.Double2String(xm,checkTOL,prec,fixW) << ' ';
      if (type==2) {
        sout() << GlobalVar.Double2String(xn,checkTOL,prec,fixW);
      } else if (type == 3) {
        sout() << GlobalVar.Double2String(xn/(tdouble)N,checkTOL,prec,fixW);
      } else if (type == 4 || type==5 ) {
        sout() << GlobalVar.Double2String(xn/(N*dx),checkTOL,prec,fixW);
      }
      if (type != 5 ) {
        sout() << ' ' << GlobalVar.Double2String(dx,checkTOL,prec,fixW);
      }
      sout()  << std::endl;;
      xe+=dx;
      xm+=dx;
      xn=0;
      if (stream_empty) break;
    }
    if (type == 5 ) {
      sout() << GlobalVar.Double2String(xm,checkTOL,prec,fixW) << ' ' << GlobalVar.Double2String(ZERO,checkTOL,prec,fixW) << std::endl;
    }
  } else {
    std::ostringstream ssV;
    ssV << "Unknown type of plot (" << type << ").";
    throw FlxException_NeglectInInteractive("FlxObjSmpPlot::task_2", ssV.str() );
  }
  isv->reset_stream();
}

FunBase* FunReadFunSmpCDF::read(bool errSerious)
{
  FlxString* flxStr = new FlxString(false,false);
  FunBase* fb = NULL;
  try {
    reader->getChar(',',false);
    fb = FunctionList->read(errSerious);
    bool inverse = false;
    if (reader->whatIsNextChar()==',') {
      reader->getChar(',');
      const std::string st = reader->getWord(true,false);
      if (st=="yes") inverse=true;
      else if (st=="no") inverse=false;
      else {
        std::ostringstream ssV;
        ssV << "Unknown keyword '" << st << "'.";
        throw FlxException_NeglectInInteractive("FunReadFunSmpCDF::read", ssV.str() );
      }
    }
    return new FunSmpCDF(flxStr,fb,inverse);
  } catch (FlxException &e) {
    FLXMSG("FunReadFunSmpCDF::read",1);
    delete flxStr;
    if (fb) delete fb;
    throw;
  }
}

FunSmpCDF::~FunSmpCDF()
{
  delete isname;
  delete val;
}

tdouble FunSmpCDF::inv_cdf(const tdouble thr, const tdouble*const tp, const tuint N)
{
  const tdouble ph = (ONE)/(2*N);
  if (thr<ZERO||thr>ONE) {
    std::ostringstream ssV;
    ssV << "Value '" << GlobalVar.Double2String(thr) << "' is not within the valid bounds [0;1]";
    throw FlxException_NeglectInInteractive("FunSmpCDF::calc_2", ssV.str() );
  }
  if (thr <= ph) return tp[0];
  if (thr >= (ONE-ph)) return tp[N-1];
  const tuint id = (thr-ph)/(2*ph);
  const tdouble pd = thr-(id+ONE/2)/N;
  #if FLX_DEBUG
    if (pd<ZERO||pd>2*ph) throw FlxException_Crude("FunSmpCDF::calc_3");
  #endif
  return tp[id]+(tp[id+1]-tp[id])*(pd/(2*ph));
}

const tdouble FunSmpCDF::calc()
{
  const std::string isname_str = isname->eval_word(true);
  FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(isname_str)));
  if (isv==NULL) {
    std::ostringstream ssV;
    ssV << "The input stream'" << isname_str << "' is not a vector-input stream.";
    throw FlxException_NeglectInInteractive("FunSmpCDF::calc_1", ssV.str() );
  }
  isv->reset_stream();
  const tulong N = isv->get_total_size();
  const tdouble thr = val->calc();
  const tdouble* const tp = isv->get_tmpPtr();
  if (inverse==false) {                // CDF
    if (thr<tp[0]) return ZERO;
    else if (thr>tp[N-1]) return ONE;
    const tuint id = flx_find_pos2(tp,N,thr);  // tp[id]<=thr
    if (id==N) return (N-ONE/2)/N;
    const tdouble pl = (id+ONE/2)/N;
    return pl+(thr-tp[id])/(tp[id+1]-tp[id])*(ONE/N);
  } else {                        // inverse CDF
    return inv_cdf(thr,tp,N);
  }
}

const std::string FunSmpCDF::write()
{
  std::ostringstream ssV;
  ssV << "cdf_smp(" << isname->write() << "," << val->write();
  if (inverse) ssV << ",yes";
  ssV << ")";
  return ssV.str();
}

FlxObjQQplot::FlxObjQQplot(const bool dolog, const std::string ostreamV, FlxString* isname, RBRV_entry_RV_base* rep)
: FlxObjOutputBase(dolog,ostreamV), isname(isname), rep(rep)
{
  
}

FlxObjQQplot::~FlxObjQQplot()
{
  delete isname;
  delete rep;
}

void FlxObjQQplot::task()
{
  rep->eval_para();
  // initializations and preparations
    const std::string isname_str = isname->eval_word(true);
    FlxIstream_vector *isv = dynamic_cast<FlxIstream_vector*>(&(data->IstreamBox.get(isname_str)));
    if (isv==NULL) {
      std::ostringstream ssV;
      ssV << "The input stream'" << isname_str << "' is not a vector-input stream.";
      throw FlxException_NeglectInInteractive("FlxObjQQplot::task_1", ssV.str() );
    }
    isv->reset_stream();
    const tulong N = isv->get_total_size();
    GlobalVar.slog(4) << "qq_plot: vector-input stream with " << GlobalVar.Double2String(isv->get_total_size()) << " entries." << std::endl;
  // the main loop
    tdouble d;
    for (tulong i=0;i<N;++i) {
      isv->get_value(d,true);
      const tdouble y = (i<N/2)?(rv_InvPhi((tdouble(i)+ONE/2)/N)):(-rv_InvPhi((N-(tdouble(i)+ONE/2))/N));
      const tdouble x = rep->transform_y2x(y);
      sout() << GlobalVar.Double2String(d) << ' ' << GlobalVar.Double2String(x) << std::endl; 
    }
  isv->reset_stream();
}


FlxObjBase* FlxObjReadQQplot::read()
{
  FlxString* flxStr = NULL;
  RBRV_entry_RV_base* rep = NULL;
  try {
    reader->getChar('(',false);
    // get the input stream
      flxStr = new FlxString(false,false);
      reader->getChar(',',false);
    // get the distribution information
      rep = RBRV_entry_read_base::read_gen_entry(false);
    reader->getChar(')',false);
    read_optionalPara(false);
    return new FlxObjQQplot(get_doLog(),get_stream(),flxStr,rep);
  } catch (FlxException &e) {
    FLXMSG("FlxObjReadQQplot::read",1);
    if (flxStr) delete flxStr;
    if (rep) delete rep;
    throw;
  }
}


FlxObjSensi_s1o_new::FlxObjSensi_s1o_new ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* Nlearn, FlxFunction* x_dim)
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), Nlearn(Nlearn), x_dim(x_dim)
{

}

void FlxObjSensi_s1o_new::task()
{
  const std::string name = nameID->eval_word(true);
  const size_t N = Nlearn->cast2tuint(false);
  const tuint M = x_dim->cast2tuint(false);
  flx_sensi_s1o* sa = new flx_sensi_s1o(name,N,M);
  sensi_s1o_box.insert(name,sa);
}

FlxObjReadSensi_s1o_new::FlxObjReadSensi_s1o_new() {
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*10000,"sensi::nlearn"));
    ParaBox.insert("nlearn", "sensi::nlearn" );
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE,"sensi::xdim"));
    ParaBox.insert("xdim", "sensi::xdim" );
}

FlxObjBase* FlxObjReadSensi_s1o_new::read()
{
  FlxString* nameID = new FlxString(false,false);
  try {
    read_optionalPara(false);
    return new FlxObjSensi_s1o_new(get_doLog(),get_stream(),nameID,get_optPara_FlxFunction("nlearn"),get_optPara_FlxFunction("xdim"));
  } catch (FlxException& e) {
    delete nameID;
    throw;
  }
}


FlxObjSensi_s1o_add::FlxObjSensi_s1o_add ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* value_x, FlxMtxConstFun* xvec, FlxFunction* value_y )
: FlxObjOutputBase(dolog,ostreamV), nameID(nameID), value_x(value_x), value_y(value_y), xvec(xvec)
{

}

void FlxObjSensi_s1o_add::task()
{
  // retrieve sensitivity object
    const std::string name = nameID->eval_word(true);
    flx_sensi_s1o& sa = sensi_s1o_box.get(name);
  // evalute x
    flxVec xvals(sa.get_x_dim());
    xvals.set_nan();
    if (xvec) {
      FlxSMtx* xvec_s = data->ConstMtxBox.get(xvec->eval(),true);
      tuint m = xvec_s->get_Ncoeff();
      if (m>sa.get_x_dim()) m = sa.get_x_dim();
      for (tuint i=0;i<m;++i) {
        xvals[i] = xvec_s->operator()(i);
      }
    } else {
      xvals[0] = value_x->calc();
    }
  // evaluate y
    const tdouble yval = value_y->calc();
  // record value
    sa.record_value(xvals,yval);
}

FlxObjSensi_s1o_add::~FlxObjSensi_s1o_add()
{
  delete nameID;
  if (value_x) {
    delete value_x;
  }
  delete value_y;
  if (xvec) {
    delete xvec;
  }
}

FlxObjBase* FlxObjReadSensi_s1o_add::read()
{
  FlxString* nameID = new FlxString(false,false);
  FlxFunction* value_x = NULL;
  FlxMtxConstFun* xvec = NULL;
  FlxFunction* value_y = NULL;
  try {
    reader->getChar('+');
    reader->getChar('=');
    reader->getChar('(');
    if (reader->whatIsNextChar()=='{') {
        reader->getChar('{',false);
        xvec = new FlxMtxConstFun(true);
        reader->getChar('}',false);
    } else {
      value_x = new FlxFunction(funReader,false);
    }
    reader->getChar(',');
    value_y = new FlxFunction(funReader,false);
    reader->getChar(')');
    read_optionalPara(false);
    return new FlxObjSensi_s1o_add(get_doLog(),get_stream(),nameID,value_x,xvec,value_y);
  } catch (FlxException& e) {
    delete nameID;
    if (value_x) delete value_x;
    if (value_y) delete value_y;
    if (xvec) delete xvec;
    throw;
  }
}

const tdouble FunSensi_s1o_eval::calc()
{
  const std::string name = nameID->eval_word(true);
  flx_sensi_s1o& sa = sensi_s1o_box.get(name);
  return sa.eval();
}

const std::string FunSensi_s1o_eval::write()
{
  std::ostringstream ssV;
  ssV << "sensi_s1o_eval(" << nameID->write() << ")";
  return ssV.str();
}

FunBase* FunReadFunSensi_s1o_eval::read( bool errSerious )
{
    FlxString* flxStr = new FlxString(false,false);
    return new FunSensi_s1o_eval(flxStr);
}



FlxObjReadSensi_s1o_dist::FlxObjReadSensi_s1o_dist() {
  AllDefParaBox->insert(new FlxOptionalParaFun(ONE*1000,"sensi::n"));
    ParaBox.insert("n", "sensi::n" );
}

FlxObjBase* FlxObjReadSensi_s1o_dist::read()
{
  // read name of vector to store
  FlxMtxConstFun* MtxConstStr = new FlxMtxConstFun(false);
  // get name of sensi Object
  FlxString* nameID = NULL;
  try {
    reader->getChar('=');
    nameID = new FlxString(false,false);
    read_optionalPara(false);
    return new FlxObjSensi_s1o_dist(get_doLog(),nameID,MtxConstStr,get_optPara_FlxFunction("n"));
  } catch (FlxException& e) {
    delete MtxConstStr;
    if (nameID) delete nameID;
    throw;
  }
}

FlxObjSensi_s1o_dist::~FlxObjSensi_s1o_dist()
{
  delete MtxConstStr;
  delete Nfun;
}

void FlxObjSensi_s1o_dist::task()
{
  // get sensitivity object
    const std::string name = nameID->eval_word(true);
    flx_sensi_s1o& sa = sensi_s1o_box.get(name);
  // get vector to store results
    const tuint N = Nfun->cast2tuint(false);
    const std::string vecName = MtxConstStr->eval();
    tdouble* vp = data->ConstMtxBox.get_Vec(N,vecName);
    flxVec svec(vp,N);
  // perform sampling procedure
    sa.eval_dist(svec,data->RndCreator);
}

