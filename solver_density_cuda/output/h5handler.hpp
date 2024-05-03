#pragma once 

//https://qiita.com/Miya42/items/600c90d09c41738535c3

#include<hdf5.h>
#include<cassert>
#include<iostream>
#include<vector>
#include<string>
#include<type_traits>

template<typename T>
constexpr bool false_v = false;

template <typename T>
hid_t get_predtype()
{
  if constexpr (std::is_same<T,char>{}){
    return H5T_NATIVE_CHAR;
  }else if constexpr (std::is_same<T,int>{}){
    return H5T_NATIVE_INT;  
  }else if constexpr (std::is_same<T,float>{}){
    return H5T_NATIVE_FLOAT;  
  }else if constexpr (std::is_same<T,double>{}){
    return H5T_NATIVE_DOUBLE;  
  }else{
    static_assert(false_v<T>,"cannot find predtype");
  }
  return 0;
}

template <typename T>
hid_t get_dspace_type()
{
  if constexpr (std::is_same<T,char>{}){
    return H5T_STD_I8LE;
  }else if constexpr (std::is_same<T,int>{}){
    return H5T_STD_I32LE;  
  }else if constexpr (std::is_same<T,float>{}){
    return H5T_IEEE_F32LE;  
  }else if constexpr (std::is_same<T,double>{}){
    return H5T_IEEE_F64LE;  
  }else{
    static_assert(false_v<T>,"cannot find dspacetype");
  }
  return 0;
}



struct h5fp
{
private:
  hid_t fid;

public:
  h5fp(std::string fname, char mode)
  { open(fname,mode); };

  ~h5fp()
  { close(); }

  void open(std::string fname, char mode)
  {
    switch (mode)
    {
    case 'w':
      fid = H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      break;
    case 'r':
      fid = H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);      
      break;    
    default:
      assert(!"[error] invalid mode from h5f::open"); 
      break;
    }
  }

  void close()
  {
    H5Fclose(fid);    
    return;
  }

  void create_group(std::string gname)
  {
    hid_t gid = H5Gcreate(fid,gname.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Gclose(gid);
    return;
  }

  template<typename T>
  void create_dataset(std::string dname,std::vector<int> shape,T data[])
  {
    size_t sz = shape.size();
    hsize_t rank = (hsize_t)sz;
    hsize_t dims[sz];
    for (size_t i = 0; i < sz; i++) dims[i] = (hsize_t)shape[i];
    hid_t dataspace = H5Screate_simple(rank,dims,NULL);
    hid_t dataset = H5Dcreate(fid,dname.c_str(),get_dspace_type<T>(),dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(dataset,get_predtype<T>(),H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return;
  }

  template<typename T>
  void update(std::string dname,std::vector<int> shape,std::vector<int> offset,T data[])
  {
    hid_t dataset = H5Dopen(fid,dname.c_str(),H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    size_t sz = shape.size();
    hsize_t rank = (hsize_t)sz;
    hsize_t dims[sz];
    hsize_t offsets[sz];
    for (size_t i = 0; i < sz; i++){ 
      dims[i] = (hsize_t)shape[i];
      offsets[i] = (hsize_t)offset[i];
    }
    hid_t dataspace_buf = H5Screate_simple(rank,dims,NULL);
    H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offsets,NULL,dims,NULL);
    H5Dwrite(dataset,get_predtype<T>(),dataspace_buf,dataspace,H5P_DEFAULT,data);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return;
  }

  template<typename T>
  void create_dataset(std::string dname,std::vector<int> shape,T fillval)
  {
    size_t sz = shape.size();
    hsize_t rank = (hsize_t)sz;
    hsize_t dims[sz];
    for (size_t i = 0; i < sz; i++) dims[i] = (hsize_t)shape[i];
    hid_t dataspace = H5Screate_simple(rank,dims,NULL);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_fill_value(plist,get_predtype<T>(),&fillval);
    hid_t dataset = H5Dcreate(fid,dname.c_str(),get_dspace_type<T>(),dataspace,H5P_DEFAULT,plist,H5P_DEFAULT);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(plist);
    return;
  }

  template<typename T>
  void read_dataset(std::string dname, T data[] )
  {
    hid_t dataset = H5Dopen(fid,dname.c_str(),H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    H5Dread(dataset, get_predtype<T>(), dataspace, dataspace, H5P_DEFAULT, data);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    return;
  }

  std::vector<hsize_t> get_shape(std::string dname)
  {
    hid_t dataset = H5Dopen(fid,dname.c_str(),H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t rank = H5Sget_simple_extent_ndims(dataspace);
    std::vector<hsize_t> dims(rank);
    H5Sget_simple_extent_dims(dataspace,dims.data(),NULL);
    H5Dclose(dataset);
    H5Sclose(dataspace);    
    return dims;
  }
};