// File       : XDMFDriver.h
// Created    : Sat Jan 25 2020 02:42:37 PM (+0100)
// Author     : Fabian Wermelinger
// Description: XDMF meta data interface
// Copyright 2020 ETH Zurich. All Rights Reserved.
#ifndef XDMFDRIVER_H_WQMJLHYN
#define XDMFDRIVER_H_WQMJLHYN

#include "Cubism/Common.h"
#include <cstdio>
#include <sstream>
#include <string>
#include <type_traits>

NAMESPACE_BEGIN(Cubism)
NAMESPACE_BEGIN(IO)

template <typename DataType, Cubism::MeshClass Class>
struct XDMFDriver {
    template <typename Mesh>
    void write(const std::string &,
               const std::string &,
               const Mesh &,
               const Cubism::EntityType,
               const size_t,
               const typename Mesh::PointType origin,
               const double) const;
};

template <typename DataType>
struct XDMFDriver<DataType, Cubism::MeshClass::Uniform> {
    template <typename Mesh>
    void write(const std::string &fname,
               const std::string &aname,
               const Mesh &mesh,
               const Cubism::EntityType entity,
               const size_t NComp,
               const typename Mesh::PointType origin,
               const double time) const
    {
        // XXX: [fabianw@mavt.ethz.ch; 2020-01-29] The ParaView XDMF reader
        // seems to be buggy for Mesh::Dim == 2 and Cubism::EntityType::Node.
        std::string topology("");
        std::string geometry("");
        if (2 == Mesh::Dim) {
            topology = "2DCoRectMesh";
            geometry = "Origin_DxDy";
        } else if (3 == Mesh::Dim) {
            topology = "3DCoRectMesh";
            geometry = "Origin_DxDyDz";
        }
        std::string data_attr("");
        if (1 == NComp) {
            data_attr = "Scalar";
        } else {
            if (Mesh::Dim == NComp) {
                data_attr = "Vector";
            } else {
                data_attr = "Tensor";
            }
        }
        std::string data_center("");
        typename Mesh::MultiIndex data_dims;
        if (entity == Cubism::EntityType::Cell) {
            data_center = "Cell";
            data_dims =
                mesh.getIndexRange(Cubism::EntityType::Cell).getExtent();
        } else if (entity == Cubism::EntityType::Node) {
            data_center = "Node";
            data_dims =
                mesh.getIndexRange(Cubism::EntityType::Node).getExtent();
        }
        std::string data_type("Float");
        if (std::is_integral<DataType>::value) {
            if (std::is_unsigned<DataType>::value) {
                data_type = "UInt";
                if (1 == sizeof(DataType)) {
                    data_type = "UChar";
                }
            } else {
                data_type = "Int";
                if (1 == sizeof(DataType)) {
                    data_type = "Char";
                }
            }
        }

        std::ostringstream mdims; // mesh dimensions
        std::ostringstream ddims; // data dimensions
        std::ostringstream orig;  // origin
        std::ostringstream spac;  // spacing
        orig.precision(16);
        spac.precision(16);
        const auto nodes =
            mesh.getIndexRange(Cubism::EntityType::Node).getExtent();
        const auto spacing = mesh.getCellSize(0);
        mdims << nodes[Mesh::Dim - 1];
        ddims << data_dims[Mesh::Dim - 1];
        orig << origin[Mesh::Dim - 1];
        spac << spacing[Mesh::Dim - 1];
        for (size_t i = 1; i < Mesh::Dim; ++i) {
            mdims << " " << nodes[Mesh::Dim - 1 - i];
            ddims << " " << data_dims[Mesh::Dim - 1 - i];
            orig << " " << origin[Mesh::Dim - 1 - i];
            spac << " " << spacing[Mesh::Dim - 1 - i];
        }
        ddims << " " << NComp;
        std::string mesh_dimZYX(mdims.str());
        std::string data_dimZYXC(ddims.str());
        std::string mesh_origin(orig.str());
        std::string mesh_spacing(spac.str());

        // XXX: [fabianw@mavt.ethz.ch; 2020-01-26] Remove path; Linux/Mac only
        const std::string basename = fname.substr(fname.find_last_of("/") + 1);

        std::FILE *xmf = 0;
        xmf = fopen((fname + ".xmf").c_str(), "w");
        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
        fprintf(xmf, "<Domain>\n");
        fprintf(xmf, "\t<Grid GridType=\"Uniform\">\n");
        fprintf(xmf, "\t\t<Time Value=\"%e\"/>\n\n", time);
        fprintf(xmf,
                "\t\t<Topology TopologyType=\"%s\" Dimensions=\"%s\"/>\n\n",
                topology.c_str(),
                mesh_dimZYX.c_str());
        fprintf(xmf, "\t\t<Geometry GeometryType=\"%s\">\n", geometry.c_str());
        fprintf(xmf,
                "\t\t\t<DataItem Name=\"Origin\" Dimensions=\"%zu\" "
                "NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n",
                Mesh::Dim);
        fprintf(xmf, "\t\t\t\t%s\n", mesh_origin.c_str());
        fprintf(xmf, "\t\t\t</DataItem>\n");
        fprintf(xmf,
                "\t\t\t<DataItem Name=\"Spacing\" Dimensions=\"%zu\" "
                "NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n",
                Mesh::Dim);
        fprintf(xmf, "\t\t\t\t%s\n", mesh_spacing.c_str());
        fprintf(xmf, "\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t</Geometry>\n\n");
        fprintf(xmf,
                "\t\t<Attribute Name=\"%s\" AttributeType=\"%s\" "
                "Center=\"%s\">\n",
                aname.c_str(),
                data_attr.c_str(),
                data_center.c_str());
        fprintf(xmf,
                "\t\t\t<DataItem Dimensions=\"%s\" NumberType=\"%s\" "
                "Precision=\"%zu\" Format=\"HDF\">\n",
                data_dimZYXC.c_str(),
                data_type.c_str(),
                sizeof(DataType));
        fprintf(xmf, "\t\t\t\t./%s:/data\n", (basename + ".h5").c_str());
        fprintf(xmf, "\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t</Attribute>\n");
        fprintf(xmf, "\t</Grid>\n");
        fprintf(xmf, "</Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
    }
};

// template <DataType>
// struct XDMFDriver<DataType, Cubism::MeshClass::Stretched> {
//     template <typename Mesh>
//     void write(const std::string &fname,
//                const Mesh &mesh,
//                const Cubism::EntityType entity,
//                const typename Mesh::PointType origin,
//                const double time) const
//     {
// TODO: [fabianw@mavt.ethz.ch; 2020-01-25]
//     }
// };

NAMESPACE_END(IO)
NAMESPACE_END(Cubism)

#endif /* XDMFDRIVER_H_WQMJLHYN */
