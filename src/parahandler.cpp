
#include "parahandler.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

ParaHandler::ParaHandler(char *model_file) {
  this->read_model_info(model_file);
}

ParaHandler::~ParaHandler() {}

void ParaHandler::skip_comments(std::istream &in,
                                std::vector<std::string> &para_str_vec,
                                std::string comment_str) {
  std::string line;
  while (std::getline(in, line)) {
    for (char &c : line) // loop string in C++11 grammar
    {
      if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n') {
        c = ' ';
      }
    }

    line.erase(0, line.find_first_not_of(
                      " ")); // delete the space at the beginning of the line
    line.erase(line.find_last_not_of(" ") +
               1); // delete the space at the ending of the line

    int n_comment_start = line.find_first_of(comment_str);
    if (n_comment_start != std::string::npos) {
      line.erase(n_comment_start); // delete comments
    }

    if (line.empty())
      continue;

    para_str_vec.push_back(line);
  }
}

void ParaHandler::read_model_info(char *model_file) {
  std::ifstream in_stream(model_file);
  assert(in_stream.good());
  int para_lines = 12; // number of parameters' lines
  std::vector<std::string> para_str_vec;
  skip_comments(in_stream, para_str_vec);
  assert(para_str_vec.size() == para_lines);
  in_stream.close();
  std::stringstream ss;

  // survey mode file
  ss << para_str_vec[0];
  ss >> survey_input;
  ss.clear();
  assert(survey_input == "survey.input");

  // Model parameters
  ss << para_str_vec[1];
  ss >> model_parameters_file;
  ss.clear();
  assert(model_parameters_file == "sigma.file");

  ss << para_str_vec[2];
  ss >> maxit;
  ss.clear();
  assert(maxit >= 1);

  ss << para_str_vec[3];
  ss >> amr_type;
  ss.clear();
  assert(amr_type == 0 || amr_type == 1);

  ss << para_str_vec[4];
  ss >> beta;
  ss.clear();
  assert(beta >= 0.0 && beta <= 1);

  ss << para_str_vec[5];
  ss >> max_dofs;
  ss.clear();
  assert(max_dofs >= 10);

  ss << para_str_vec[6];
  ss >> linear_solver;
  ss.clear();
  // assert(linear_solver == "mumps" || linear_solver == "ilu");

  // Mesh parameters
  ss << para_str_vec[7];
  ss >> save_amr_mesh;
  ss.clear();

  ss << para_str_vec[8];
  ss >> mesh_file;
  ss.clear();

  ss << para_str_vec[9];
  ss >> if_compute_P2;
  ss.clear();
  assert(if_compute_P2 == "true" || if_compute_P2 == "false");

  ss << para_str_vec[10];
  ss >> if_compute_J;
  ss.clear();
  assert(if_compute_J == "true" || if_compute_J == "false");

  ss << para_str_vec[11];
  ss >> if_use_improved_AMR;
  ss.clear();
  assert(if_use_improved_AMR == "true" || if_use_improved_AMR == "false");

  // sorted_survey = "sorted_" + survey_input;
  unsigned int n_data = 0;
  int mode = 0;

  std::ifstream survey_in_stream(survey_input);
  assert(survey_in_stream.good());
  survey_in_stream >> n_data // nunmber of measurments
      >> mode;               // measurment mode
  Vertex last_ts(-1, -1, -1);
  Vertex ts(-1, -1, -1);
  Vertex tm(-1, -1, -1);
  Vertex ts2(-1, -1, -1);
  Vertex ts3(-1, -1, -1);
  Vertex ts4(-1, -1, -1);
  std::vector<Vertex> temp_sites;
  for (int i = 0; i < n_data; i++) {
    if (mode == 11) {
      survey_in_stream >> ts(0) >> ts(1) >> ts(2);
      survey_in_stream >> tm(0) >> tm(1) >> tm(2);
      if (i == 0) {
        last_ts = ts;
        sources.push_back(ts);
        temp_sites.push_back(tm);
      } else if (i == n_data - 1) {
        temp_sites.push_back(tm);
        sites.push_back(temp_sites);
        temp_sites.clear();
      } else {
        if (ts(0) == last_ts(0) && ts(1) == last_ts(1) && ts(2) == last_ts(2)) {
          temp_sites.push_back(tm);
        } else {
          last_ts = ts;
          sources.push_back(ts);
          sites.push_back(temp_sites);
          temp_sites.clear();
          temp_sites.push_back(tm);
        }
      }
    }

    else if (mode == 41) {
      // AB-CD situation
      survey_in_stream >> ts(0) >> ts(1) >> ts(2) >> ts2(0) >> ts2(1) >>
          ts2(2) >> ts3(0) >> ts3(1) >> ts3(2) >> ts4(0) >> ts4(1) >> ts4(2);
      survey_in_stream >> tm(0) >> tm(1) >> tm(2);
      if (i == 0) {
        sources.push_back(ts);
        sources.push_back(ts2);
        sources.push_back(ts3);
        sources.push_back(ts4);
        temp_sites.push_back(tm);
      } else if (i == n_data - 1) {
        temp_sites.push_back(tm);
        sites.push_back(temp_sites);
        sites.push_back(temp_sites);
        sites.push_back(temp_sites);
        sites.push_back(temp_sites);
      } else {
        temp_sites.push_back(tm);
      }

    } else {
      std::cout << "Unsopported mode temporarily!\n";
      std::abort();
    }
  }
  survey_in_stream.close();
  source_number = sources.size();
  assert(source_number > 0);

  // Read conductivity model
  std::ifstream cond_in_stream(model_parameters_file.c_str());
  assert(cond_in_stream.good());
  cond_in_stream >> n_regions;
  for (int i = 0; i < n_regions; i++) {
    int marker;
    Vector main_cond(3);
    main_cond = 0.0;
    Vector angle(3);
    angle = 0.0;
    // double permeability;
    cond_in_stream >> marker >> main_cond[0] >> main_cond[1] >> main_cond[2] >>
        angle[0] >> angle[1] >> angle[2]; // >> permeability;
    marker_vec.push_back(marker);
    region2conductivity[marker] = cal_conductivity(main_cond, angle);
  }
  cond_in_stream.close();
  // check mesh file stream
  std::ifstream msh_stream(mesh_file);
  assert(msh_stream.good());
  msh_stream.close();
}

DenseMatrix ParaHandler::cal_conductivity(Vector &main_cond, Vector &_angle) {
  Vector angle(3);
  angle = 0.0;
  angle[0] = _angle[0] * (DC::PI / 180);
  angle[1] = _angle[1] * (DC::PI / 180);
  angle[2] = _angle[2] * (DC::PI / 180);
  DenseMatrix sigma(3);
  sigma = 0.0;
  sigma(0, 0) = main_cond[0];
  sigma(1, 1) = main_cond[1];
  sigma(2, 2) = main_cond[2];

  DenseMatrix R_1(3), R_2(3), R_3(3), R(3);
  R_1 = 0.0;
  R_2 = 0.0;
  R_3 = 0.0;
  R_1(0, 0) = cos(angle[0]);
  R_1(0, 1) = -sin(angle[0]);
  R_1(0, 2) = 0;
  R_1(1, 0) = sin(angle[0]);
  R_1(1, 1) = cos(angle[0]);
  R_1(1, 2) = 0;
  R_1(2, 0) = 0;
  R_1(2, 1) = 0;
  R_1(2, 2) = 1;
  R_2(0, 0) = 1;
  R_2(0, 1) = 0;
  R_2(0, 2) = 0;
  R_2(1, 0) = 0;
  R_2(1, 1) = cos(angle[1]);
  R_2(1, 2) = -sin(angle[1]);
  R_2(2, 0) = 0;
  R_2(2, 1) = sin(angle[1]);
  R_2(2, 2) = cos(angle[1]);
  R_3(0, 0) = cos(angle[2]);
  R_3(0, 1) = -sin(angle[2]);
  R_3(0, 2) = 0;
  R_3(1, 0) = sin(angle[2]);
  R_3(1, 1) = cos(angle[2]);
  R_3(1, 2) = 0;
  R_3(2, 0) = 0;
  R_3(2, 1) = 0;
  R_3(2, 2) = 1;

  DenseMatrix temp1(3);
  temp1 = 0.0;
  Mult_DenseMatrix3(R_1, R_2, temp1);
  Mult_DenseMatrix3(temp1, R_3, R);

  DenseMatrix R_sigma_RT(3), temp2(3), temp3(3);
  R_sigma_RT = 0.0;
  temp2 = 0.0;
  temp3 = R;
  temp3.Transpose();
  Mult_DenseMatrix3(R, sigma, temp2);
  Mult_DenseMatrix3(temp2, temp3, R_sigma_RT);
  return R_sigma_RT;
}

DenseMatrix ParaHandler::get_elem_conductivity(int marker) {
  assert(!region2conductivity.empty());
  std::map<int, DenseMatrix>::iterator it = region2conductivity.find(marker);
  if (it == region2conductivity.end()) {
    std::cout << "not found marker: " << marker << std::endl;
  }
  assert(it != region2conductivity.end());
  return (*it).second;
}

Vector ParaHandler::get_sources_center() {
  Vector gpoint(3);
  gpoint = 0.0;
  for (int i = 0; i < sources.size(); i++) {
    for (int j = 0; j < 3; j++) {
      gpoint[j] += sources[i](j);
    }
  }
  gpoint /= sources.size();
  return gpoint;
}

void ParaHandler::preProcess() {
  for (int i = 0; i < sites.size(); i++) {
    for (int j = 0; j < sites[i].size(); j++) {
      s_plus_m.push_back(sites[i][j]);
    }
  }
  remove_duplicate_vertex(s_plus_m);
}

void ParaHandler::find_point_tets(ParMesh *pmesh, std::vector<Vertex> &points,
                                  Array<int> &find_tets,
                                  std::string find_points_by) {
  // Find the potential tets which contain sites
  int npts = points.size();
  DenseMatrix point_mat(3, npts);
  for (int j = 0; j < npts; j++) {
    point_mat(0, j) = points[j](0);
    point_mat(1, j) = points[j](1);
    point_mat(2, j) = points[j](2);
  }
  Array<int> temp_tet_id;
  Array<IntegrationPoint> ips;
  pmesh->FindPoints(point_mat, temp_tet_id, ips);
  for (int j = 0; j < temp_tet_id.Size(); j++) {
    int tid = temp_tet_id[j];
    if (tid == -1) {
      std::cout << "measurements point not found! please check your parameters";
    } else {
      find_tets.Append(tid);
    }
  }
  assert(find_tets.Size() == npts);
}
