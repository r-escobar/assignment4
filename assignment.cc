#include <dirent.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>
#include <jpeglib.h>

// interleaved RGB image struct RGB RGB RGB, row major:
// RGBRGBRGB
// RGBRGBRGB
// RGBRGBRGB
// above: example 3 x 3 image.
// 8 bits per channel.
struct Image {
  unsigned char* bytes;
  int width;
  int height;
};

// new code
struct Bone
{
  int id;
  int parent;
  float dx;
  float dy;
  float dz;
  float length;
  glm::vec3 tangent;
  glm::vec3 normal;
  glm::vec3 binormal;
  glm::mat4 translation;
  glm::mat4 rotation;
  glm::vec4 origin;
  glm::vec4 endpoint;
};

std::vector<Bone> bone_vector;

std::vector<glm::vec4> ogre_vertices;
std::vector<glm::uvec3> ogre_faces;

std::vector<glm::vec4> bone_vertices;
std::vector<glm::uvec2> bone_lines;
// end new code

int window_width = 800, window_height = 600;
const std::string window_title = "Virtual Mannequin";

const float kNear = 0.0001f;
const float kFar = 1000.0f;
const float kFov = 45.0f;
float aspect = static_cast<float>(window_width) / window_height;

// Floor info.
const float eps = 0.5 * (0.025 + 0.0175);
const float kFloorXMin = -100.0f;
const float kFloorXMax = 100.0f;
const float kFloorZMin = -100.0f;
const float kFloorZMax = 100.0f;
const float kFloorY = -0.75617 - eps;

const float kCylinderRadius = 0.05f;

enum {
  kMouseModeCamera,
  kMouseModeSkeleton,
  kNumMouseModes
};
int current_mouse_mode = 0;

// VBO and VAO descriptors.

// We have these VBOs available for each VAO.
enum {
  kVertexBuffer,
  kIndexBuffer,
  kNumVbos
};

// These are our VAOs.
enum {
  kFloorVao,
  kOgreVao,
  kSkeletonVao,
  kNumVaos
};

GLuint array_objects[kNumVaos];  // This will store the VAO descriptors.
GLuint buffer_objects[kNumVaos][kNumVbos];  // These will store VBO descriptors.

float last_x = 0.0f, last_y = 0.0f, current_x = 0.0f, current_y = 0.0f;
bool drag_state = false;
int current_button = -1;
float camera_distance = 0.7;
float pan_speed = 0.1f;
float roll_speed = 0.1f;
float rotation_speed = 0.05f;
float zoom_speed = 0.1f;
glm::vec3 eye = glm::vec3(0.0f, 0.0f, camera_distance);
glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
glm::vec3 look = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 tangent = glm::cross(up, look);
glm::vec3 center = eye + camera_distance * look;
glm::mat3 orientation = glm::mat3(tangent, up, look);
bool fps_mode = false;

glm::mat4 view_matrix = glm::lookAt(eye, center, up);
glm::mat4 projection_matrix =
    glm::perspective((float)(kFov * (M_PI / 180.0f)), aspect, kNear, kFar);
glm::mat4 model_matrix = glm::mat4(1.0f);
glm::mat4 floor_model_matrix = glm::mat4(1.0f);

const char* vertex_shader =
    "#version 330 core\n"
    "uniform vec4 light_position;"
    "in vec4 vertex_position;"
    "out vec4 vs_light_direction;"
    "void main() {"
    "gl_Position = vertex_position;"
    "vs_light_direction = light_position - gl_Position;"
    "}";

const char* geometry_shader =
    "#version 330 core\n"
    "layout (triangles) in;"
    "layout (triangle_strip, max_vertices = 3) out;"
    "uniform mat4 projection;"
    "uniform mat4 model;"
    "uniform mat4 view;"
    "uniform vec4 light_position;"
    "in vec4 vs_light_direction[];"
    "out vec4 face_normal;"
    "out vec4 light_direction;"
    "out vec4 world_position;"
    "void main() {"
    "int n = 0;"
    "vec3 a = gl_in[0].gl_Position.xyz;"
    "vec3 b = gl_in[1].gl_Position.xyz;"
    "vec3 c = gl_in[2].gl_Position.xyz;"
    "vec3 u = normalize(b - a);"
    "vec3 v = normalize(c - a);"
    "face_normal = normalize(vec4(normalize(cross(u, v)), 0.0));"
    "for (n = 0; n < gl_in.length(); n++) {"
    "light_direction = normalize(vs_light_direction[n]);"
    "world_position = gl_in[n].gl_Position;"
    "gl_Position = projection * view * model * gl_in[n].gl_Position;"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";

const char* floor_fragment_shader =
    "#version 330 core\n"
    "in vec4 face_normal;"
    "in vec4 light_direction;"
    "in vec4 world_position;"
    "out vec4 fragment_color;"
    "void main() {"
    "vec4 pos = world_position;"
    "float check_width = 0.25;"
    "float i = floor(pos.x / check_width);"
    "float j  = floor(pos.z / check_width);"
    "vec3 color = mod(i + j, 2) * vec3(1.0, 1.0, 1.0);"
    "float dot_nl = dot(normalize(light_direction), normalize(face_normal));"
    "dot_nl = clamp(dot_nl, 0.0, 1.0);"
    "color = clamp(dot_nl * color, 0.0, 1.0);"
    "fragment_color = vec4(color, 1.0);"
    "}";


const char* skeleton_vertex_shader =
    "#version 330 core\n"
    "in vec4 vertex_position;"
    "void main() {"
    "gl_Position = vertex_position;"
    "}";

const char* skeleton_geometry_shader =
    "#version 330 core\n"
    "layout (lines) in;"
    "layout (line_strip, max_vertices = 2) out;"
    "uniform mat4 projection;"
    "uniform mat4 model;"
    "uniform mat4 view;"
    "out vec4 world_position;"
    "void main() {"
    "int n = 0;"
    "vec3 a = gl_in[0].gl_Position.xyz;"
    "vec3 b = gl_in[1].gl_Position.xyz;"
    "for (n = 0; n < gl_in.length(); n++) {"
    "world_position = gl_in[n].gl_Position;"
    "gl_Position = projection * view * model * gl_in[n].gl_Position;"
    "EmitVertex();"
    "}"
    "EndPrimitive();"
    "}";

const char* skeleton_fragment_shader =
    "#version 330 core\n"
    "in vec4 face_normal;"
    "in vec4 light_direction;"
    "in vec4 world_position;"
    "out vec4 fragment_color;"
    "void main() {"
    "vec4 pos = world_position;"
    "fragment_color = vec4(0.0, 0.0, 1.0, 1.0);"
    "}";


const char* OpenGlErrorToString(GLenum error) {
  switch (error) {
    case GL_NO_ERROR:
      return "GL_NO_ERROR";
      break;
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM";
      break;
    case GL_INVALID_VALUE:
      return "GL_INVALID_VALUE";
      break;
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION";
      break;
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY";
      break;
    default:
      return "Unknown Error";
      break;
  }
  return "Unicorns Exist";
}

#define CHECK_SUCCESS(x) \
  if (!(x)) {            \
    glfwTerminate();     \
    exit(EXIT_FAILURE);  \
  }

#define CHECK_GL_SHADER_ERROR(id)                                           \
  {                                                                         \
    GLint status = 0;                                                       \
    GLint length = 0;                                                       \
    glGetShaderiv(id, GL_COMPILE_STATUS, &status);                          \
    glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);                         \
    if (!status) {                                                          \
      std::string log(length, 0);                                           \
      glGetShaderInfoLog(id, length, nullptr, &log[0]);                     \
      std::cerr << "Line :" << __LINE__ << " OpenGL Shader Error: Log = \n" \
                << &log[0];                                                 \
      glfwTerminate();                                                      \
      exit(EXIT_FAILURE);                                                   \
    }                                                                       \
  }

#define CHECK_GL_PROGRAM_ERROR(id)                                           \
  {                                                                          \
    GLint status = 0;                                                        \
    GLint length = 0;                                                        \
    glGetProgramiv(id, GL_LINK_STATUS, &status);                             \
    glGetProgramiv(id, GL_INFO_LOG_LENGTH, &length);                         \
    if (!status) {                                                           \
      std::string log(length, 0);                                            \
      glGetProgramInfoLog(id, length, nullptr, &log[0]);                     \
      std::cerr << "Line :" << __LINE__ << " OpenGL Program Error: Log = \n" \
                << &log[0];                                                  \
      glfwTerminate();                                                       \
      exit(EXIT_FAILURE);                                                    \
    }                                                                        \
  }

#define CHECK_GL_ERROR(statement)                                             \
  {                                                                           \
    { statement; }                                                            \
    GLenum error = GL_NO_ERROR;                                               \
    if ((error = glGetError()) != GL_NO_ERROR) {                              \
      std::cerr << "Line :" << __LINE__ << " OpenGL Error: code  = " << error \
                << " description =  " << OpenGlErrorToString(error);          \
      glfwTerminate();                                                        \
      exit(EXIT_FAILURE);                                                     \
    }                                                                         \
  }

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  size_t count = std::min(v.size(), static_cast<size_t>(10));
  for (size_t i = 0; i < count; ++i) os << i << " " << v[i] << "\n";
  os << "size = " << v.size() << "\n";
  return os;
}

namespace glm {
std::ostream& operator<<(std::ostream& os, const glm::vec2& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::vec4& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::mat4& v) {
  os << glm::to_string(v);
  return os;
}

std::ostream& operator<<(std::ostream& os, const glm::mat3& v) {
  os << glm::to_string(v);
  return os;
}
}  // namespace glm

void LoadObj(const std::string& file, std::vector<glm::vec4>& vertices,
             std::vector<glm::uvec3>& indices) {
  std::ifstream in(file);
  int i = 0, j = 0;
  glm::vec4 vertex = glm::vec4(0.0, 0.0, 0.0, 1.0);
  glm::uvec3 face_indices = glm::uvec3(0, 0, 0);
  while (in.good()) {
    char c = in.get();
    switch (c) {
      case 'v':
        in >> vertex[0] >> vertex[1] >> vertex[2];
        vertices.push_back(vertex);
        break;
      case 'f':
        in >> face_indices[0] >> face_indices[1] >> face_indices[2];
        face_indices -= 1;
        indices.push_back(face_indices);
        break;
      default:
        break;
    }
  }
  in.close();
}

glm::mat4 coordMatrix(Bone* currBone, bool origin) {
  if(currBone->parent == 0) {
    if(origin)
      return currBone->translation;
    else
      return currBone->translation * currBone->rotation;
  } else {
    if(origin)
      return coordMatrix(&bone_vector[currBone->parent], false) * currBone->translation;
    else
     return coordMatrix(&bone_vector[currBone->parent], false) * currBone->translation * currBone->rotation;
  }
}

void calculateEndpoints(Bone* currBone) {
  
  if(currBone->parent == -1) {
    currBone->origin = glm::vec4(currBone->dx, currBone->dy, currBone->dz, 1.000000);
    currBone->endpoint = currBone->origin;
  } else {
    std::cout << coordMatrix(currBone, true) << "\n\n";
    currBone->origin = coordMatrix(currBone, true) * glm::vec4(0.000000, 0.000000, 0.000000, 1.000000);
    currBone->endpoint = coordMatrix(currBone, false) * glm::vec4(currBone->length, 0.000000, 0.000000, 1.000000);
  }

  //if(currBone->parent == 0) {
    int firstIndex = bone_vertices.size();
    bone_vertices.push_back(currBone->origin);
    int secondIndex = bone_vertices.size();
    bone_vertices.push_back(currBone->endpoint);  

    bone_lines.push_back(glm::uvec2(firstIndex, secondIndex));
  //}
}

glm::mat4 calculateRotation(Bone* currBone) {
    glm::mat4 tnb = glm::mat4(glm::vec4(currBone->tangent.x, currBone->tangent.y, currBone->tangent.z, 0.000000),
               glm::vec4(currBone->normal.x, currBone->normal.y, currBone->normal.z, 0.000000),
               glm::vec4(currBone->binormal.x, currBone->binormal.y, currBone->binormal.z, 0.000000),
               glm::vec4(0.000000, 0.000000, 0.000000, 1.000000));
    return glm::inverse(tnb);
}

void LoadBones(const std::string& file)
{

    // parse ogre-skeleton.bf file into bone structs
    std::ifstream in(file);
    int id, parent;
    float x, y, z;
    while (in >> id >> parent >> x >> y >> z) 
    {
      
      Bone newBone;
      newBone.id = id;
      std::cout << "newBone id is " << newBone.id << "\n";
      newBone.parent = parent;
      std::cout << "newBone parent is " << newBone.parent << "\n";
      newBone.dx = x;
      std::cout << "newBone dx is " << newBone.dx << "\n";
      newBone.dy = y;
      std::cout << "newBone dy is " << newBone.dy << "\n";
      newBone.dz = z;
      std::cout << "newBone dz is " << newBone.dz << "\n";
      bone_vector.push_back(newBone);

    }

    for(int i = 0; i < bone_vector.size(); i++)
    {
      if(bone_vector[i].parent > -1)
      {
        glm::vec3 bone_point = glm::vec3(bone_vector[i].dx, bone_vector[i].dy, bone_vector[i].dz);
        glm::vec3 parent_point = glm::vec3(bone_vector[bone_vector[i].parent].dx, bone_vector[bone_vector[i].parent].dy, bone_vector[bone_vector[i].parent].dz);
        bone_vector[i].tangent = glm::vec3(bone_point);
        bone_vector[i].length = glm::length(bone_vector[i].tangent);
        bone_vector[i].tangent = glm::normalize(bone_vector[i].tangent);

        bone_vector[i].translation = glm::mat4(glm::vec4(1.000000, 0.000000, 0.000000, 0.000000),
               glm::vec4(0.000000, 1.000000, 0.000000, 0.000000),
               glm::vec4(0.000000, 0.000000, 1.000000, 0.000000),
               glm::vec4(bone_vector[bone_vector[i].parent].length, 0.000000, 0.000000, 1.000000));



        glm::vec3 v = bone_vector[i].tangent;


        std::vector<float> offset_coords;
        offset_coords.push_back(v.x);
        offset_coords.push_back(v.y);
        offset_coords.push_back(v.z);
        std::vector<float>::iterator result =  std::min_element(std::begin(offset_coords), std::end(offset_coords));
        float mincoord = *result;

        if(mincoord == offset_coords[0]) {
          v = glm::vec3(1.000000, 0.000000, 0.000000);
        } else if(mincoord == offset_coords[1]) {
          v = glm::vec3(0.000000, 1.000000, 0.000000);
        } else if(mincoord == offset_coords[2]){
          v = glm::vec3(0.000000, 0.000000, 1.000000);
        }
        glm::vec3 normal = glm::normalize(glm::cross(bone_vector[i].tangent, v));
        glm::vec3 binormal = glm::normalize(glm::cross(bone_vector[i].tangent, normal));

        bone_vector[i].rotation = glm::mat4(glm::vec4(bone_vector[i].tangent.x, bone_vector[i].tangent.y, bone_vector[i].tangent.z, 0.000000),
             glm::vec4(normal.x, normal.y, normal.z, 0.000000),
             glm::vec4(binormal.x, binormal.y, binormal.z, 0.000000),
             glm::vec4(0.000000, 0.000000, 0.000000, 1.000000));


        bone_vector[i].normal = normal;
        bone_vector[i].binormal = binormal;


      } else if (bone_vector[i].parent == -1){

        bone_vector[i].translation = glm::mat4(glm::vec4(1.000000, 0.000000, 0.000000, 0.000000),
               glm::vec4(0.000000, 1.000000, 0.000000, 0.000000),
               glm::vec4(0.000000, 0.000000, 1.000000, 0.000000),
               glm::vec4(bone_vector[i].dx, bone_vector[i].dy, bone_vector[i].dz, 1.000000));


        bone_vector[i].rotation = glm::mat4();

        bone_vector[i].tangent = glm::vec3(1,0,0);
        bone_vector[i].normal = glm::vec3(0,1,0);
        bone_vector[i].binormal = glm::vec3(0,0,1);

      }
      // } else {
      //   glm::vec3 bone_point = glm::vec3(bone_vector[i].dx, bone_vector[i].dy, bone_vector[i].dz);
      //   glm::vec3 parent_point = glm::vec3(bone_vector[bone_vector[i].parent].dx, bone_vector[bone_vector[i].parent].dy, bone_vector[bone_vector[i].parent].dz);
      //   bone_vector[i].tangent = glm::vec3(bone_point - parent_point);
      //   bone_vector[i].length = glm::length(bone_vector[i].tangent);
      //   bone_vector[i].tangent = glm::normalize(bone_point);

        
      //   bone_vector[i].translation = bone_vector[0].translation;





      // }


       







        // // Calculate rotation matrix
        // std::vector<float> offset_coords;
        // if(bone_vector[i].parent <= 0) {
        //   offset_coords.push_back(std::abs(bone_vector[i].dx));
        //   offset_coords.push_back(std::abs(bone_vector[i].dy));
        //   offset_coords.push_back(std::abs(bone_vector[i].dz));
        // } else {
        //   offset_coords.push_back(std::abs(bone_vector[i].dx - bone_vector[bone_vector[i].parent].dx));
        //   offset_coords.push_back(std::abs(bone_vector[i].dy - bone_vector[bone_vector[i].parent].dy));
        //   offset_coords.push_back(std::abs(bone_vector[i].dz - bone_vector[bone_vector[i].parent].dz));
        // }
        // std::vector<float>::iterator result =  std::min_element(std::begin(offset_coords), std::end(offset_coords));
        // float mincoord = *result;
        // //std::cout << "Smallest offset coordinate: " << mincoord << "\n";
        // glm::vec3 offset_vector;
        // if(bone_vector[i].parent <= 0) {
        //   offset_vector = glm::normalize(glm::vec3(bone_vector[i].dx, bone_vector[i].dy, bone_vector[i].dz));      
        // } else {
        //   offset_vector = glm::normalize(glm::vec3(bone_vector[i].dx - bone_vector[bone_vector[i].parent].dx, bone_vector[i].dy - bone_vector[bone_vector[i].parent].dy, bone_vector[i].dz - bone_vector[bone_vector[i].parent].dz));
        // }
        // glm::vec3 u;
        // if(mincoord == offset_coords[0]) {
        //   u = glm::cross(offset_vector, glm::vec3(1.000000, 0.000000, 0.000000));
        // } else if(mincoord == offset_coords[1]) {
        //   u = glm::cross(offset_vector, glm::vec3(0.000000, 1.000000, 0.000000));
        // } else if(mincoord == offset_coords[2]){
        //   u = glm::cross(offset_vector, glm::vec3(0.000000, 0.000000, 1.000000));
        // } else {
        //   std::cout << "ERROR: min coordinate is not in offset vector!\n";
        // }
        // u = glm::normalize(u);
        // glm::vec3 v = glm::normalize(glm::cross(offset_vector, u));

        // bone_vector[i].rotation = glm::mat4(glm::vec4(offset_vector.x, offset_vector.y, offset_vector.z, 0.000000),
        //        glm::vec4(u.x, u.y, u.z, 0.000000),
        //        glm::vec4(v.x, v.y, v.z, 0.000000),
        //        glm::vec4(0.000000, 0.000000, 0.000000, 1.000000));

        // bone_vector[i].tangent = offset_vector;
        // bone_vector[i].normal = u;
        // bone_vector[i].binormal = v;

        if(bone_vector[i].parent >= 0) {
          bone_vector[i].rotation = calculateRotation(&bone_vector[bone_vector[i].parent]) * bone_vector[i].rotation;
        }

        //if(i < 10)
          //std::cout << "bone_vector[" << i << "] rotation is: \n" << bone_vector[i].rotation << "\n\n";

      //if(bone_vector[i].parent == 0) {
      std::cout << "bone_vector[" << i << "] length is " << bone_vector[i].length << "\n";

      std::cout << "Translation matrix: " << bone_vector[i].translation << "\n\n";
      std::cout << "Rotation matrix: " << bone_vector[i].rotation << "\n\n";    


      calculateEndpoints(&bone_vector[i]);
      std::cout << "bone_vector[" << i << "] origin is " << bone_vector[i].origin << "\n";
      std::cout << "bone_vector[" << i << "] endpoint is " << bone_vector[i].endpoint << "\n";
      std::cout << "SANITY CHECK: length between endpoint and origin is: " << glm::length(bone_vector[i].endpoint - bone_vector[i].origin) << "\n";
      std::cout << "\n";
      //}
    }

}

void SaveJPEG(const std::string& filename, int image_width, int image_height,
              const unsigned char* pixels) {
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE* outfile;
  JSAMPROW row_pointer[1];
  int row_stride;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  CHECK_SUCCESS((outfile = fopen(filename.c_str(), "wb")) != NULL)

  jpeg_stdio_dest(&cinfo, outfile);

  cinfo.image_width = image_width;
  cinfo.image_height = image_height;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, 100, true);
  jpeg_start_compress(&cinfo, true);

  row_stride = image_width * 3;

  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = const_cast<unsigned char*>(
        &pixels[(cinfo.image_height - 1 - cinfo.next_scanline) * row_stride]);
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }

  jpeg_finish_compress(&cinfo);
  fclose(outfile);

  jpeg_destroy_compress(&cinfo);
}

bool LoadJPEG(const std::string& file_name, Image* image) {
  FILE* file = fopen(file_name.c_str(), "rb");
  struct jpeg_decompress_struct info;
  struct jpeg_error_mgr err;

  info.err = jpeg_std_error(&err);
  jpeg_create_decompress(&info);

  CHECK_SUCCESS(file != NULL);

  jpeg_stdio_src(&info, file);
  jpeg_read_header(&info, true);
  jpeg_start_decompress(&info);

  image->width = info.output_width;
  image->height = info.output_height;

  int channels = info.num_components;
  long size = image->width * image->height * 3;

  image->bytes = new unsigned char[size];

  int a = (channels > 2 ? 1 : 0);
  int b = (channels > 2 ? 2 : 0);
  std::vector<unsigned char> scan_line(image->width * channels, 0);
  unsigned char* p1 = &scan_line[0];
  unsigned char** p2 = &p1;
  unsigned char* out_scan_line = &image->bytes[0];
  while (info.output_scanline < info.output_height) {
    jpeg_read_scanlines(&info, p2, 1);
    for (int i = 0; i < image->width; ++i) {
      out_scan_line[3 * i] = scan_line[channels * i];
      out_scan_line[3 * i + 1] = scan_line[channels * i + a];
      out_scan_line[3 * i + 2] = scan_line[channels * i + b];
    }
    out_scan_line += image->width * 3;
  }
  jpeg_finish_decompress(&info);
  fclose(file);
  return true;
}

void ErrorCallback(int error, const char* description) {
  std::cerr << "GLFW Error: " << description << "\n";
}

void KeyCallback(GLFWwindow* window, int key, int scancode, int action,
                 int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
  else if (key == GLFW_KEY_W && action != GLFW_RELEASE) {
    if (fps_mode)
      eye -= zoom_speed * look;
    else
      camera_distance -= zoom_speed;
  } else if (key == GLFW_KEY_S && action != GLFW_RELEASE) {
    if (fps_mode)
      eye += zoom_speed * look;
    else
      camera_distance += zoom_speed;
  } else if (key == GLFW_KEY_A && action != GLFW_RELEASE) {
    if (fps_mode)
      eye -= pan_speed * tangent;
    else
      center -= pan_speed * tangent;
  } else if (key == GLFW_KEY_D && action != GLFW_RELEASE) {
    if (fps_mode)
      eye += pan_speed * tangent;
    else
      center += pan_speed * tangent;
  } else if (key == GLFW_KEY_LEFT && action != GLFW_RELEASE) {
    if (current_mouse_mode == kMouseModeCamera) {
      glm::mat3 rotation = glm::mat3(glm::rotate(-roll_speed, look));
      orientation = rotation * orientation;
      tangent = glm::column(orientation, 0);
      up = glm::column(orientation, 1);
      look = glm::column(orientation, 2);
    } else {
      // help me
    }
  } else if (key == GLFW_KEY_RIGHT && action != GLFW_RELEASE) {
    if (current_mouse_mode == kMouseModeCamera) {
      glm::mat3 rotation = glm::mat3(glm::rotate(roll_speed, look));
      orientation = rotation * orientation;
      tangent = glm::column(orientation, 0);
      up = glm::column(orientation, 1);
      look = glm::column(orientation, 2);
    } else {
      // help me
    }
  } else if (key == GLFW_KEY_DOWN && action != GLFW_RELEASE) {
    if (fps_mode)
      eye -= pan_speed * up;
    else
      center -= pan_speed * up;
  } else if (key == GLFW_KEY_UP && action != GLFW_RELEASE) {
    if (fps_mode)
      eye += pan_speed * up;
    else
      center += pan_speed * up;
  } else if (key == GLFW_KEY_C && action != GLFW_RELEASE) {
    fps_mode = !fps_mode;
  } else if (key == GLFW_KEY_M && action != GLFW_RELEASE) {
    current_mouse_mode = (current_mouse_mode + 1) % kNumMouseModes;
  } else if (key == GLFW_KEY_J && action != GLFW_RELEASE) {
    std::vector<unsigned char> pixels(3 * window_width * window_height, 0);
    CHECK_GL_ERROR(glReadPixels(0, 0, window_width, window_height, GL_RGB,
                                GL_UNSIGNED_BYTE, &pixels[0]));
    std::string filename = "capture.jpg";
    std::cout << "Encoding and saving to file '" + filename + "'\n";
    SaveJPEG(filename, window_width, window_height, &pixels[0]);
  }
}

void MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y) {
  last_x = current_x;
  last_y = current_y;
  current_x = mouse_x;
  current_y = window_height - mouse_y;
  float delta_x = current_x - last_x;
  float delta_y = current_y - last_y;
  if (sqrt(delta_x * delta_x + delta_y * delta_y) < 1e-15) return;
  glm::vec3 mouse_direction = glm::normalize(glm::vec3(delta_x, delta_y, 0.0f));
  glm::vec2 mouse_start = glm::vec2(last_x, last_y);
  glm::vec2 mouse_end = glm::vec2(current_x, current_y);
  glm::uvec4 viewport = glm::uvec4(0, 0, window_width, window_height);
  if (drag_state && current_button == GLFW_MOUSE_BUTTON_LEFT) {
    if (current_mouse_mode == kMouseModeCamera) {
      glm::vec3 axis = glm::normalize(
          orientation * glm::vec3(mouse_direction.y, -mouse_direction.x, 0.0f));
      orientation =
          glm::mat3(glm::rotate(rotation_speed, axis) * glm::mat4(orientation));
      tangent = glm::column(orientation, 0);
      up = glm::column(orientation, 1);
      look = glm::column(orientation, 2);
    } else {
      // help me
    }
  } else {
    // maybe put some code here
  }
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
  drag_state = (action == GLFW_PRESS);
  current_button = button;
}

int main(int argc, char* argv[]) {
  if (!glfwInit()) exit(EXIT_FAILURE);
  glfwSetErrorCallback(ErrorCallback);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_SAMPLES, 4);
  GLFWwindow* window = glfwCreateWindow(window_width, window_height,
                                        &window_title[0], nullptr, nullptr);
  CHECK_SUCCESS(window != nullptr);

  glfwMakeContextCurrent(window);
  glewExperimental = GL_TRUE;
  CHECK_SUCCESS(glewInit() == GLEW_OK);
  glGetError();  // clear GLEW's error for it

  glfwSetKeyCallback(window, KeyCallback);
  glfwSetCursorPosCallback(window, MousePosCallback);
  glfwSetMouseButtonCallback(window, MouseButtonCallback);
  glfwSwapInterval(1);
  const GLubyte* renderer = glGetString(GL_RENDERER);  // get renderer string
  const GLubyte* version = glGetString(GL_VERSION);    // version as a string
  std::cout << "Renderer: " << renderer << "\n";
  std::cout << "OpenGL version supported:" << version << "\n";

  std::vector<glm::vec4> floor_vertices;
  std::vector<glm::uvec3> floor_faces;
  floor_vertices.push_back(glm::vec4(kFloorXMin, kFloorY, kFloorZMax, 1.0f));
  floor_vertices.push_back(glm::vec4(kFloorXMax, kFloorY, kFloorZMax, 1.0f));
  floor_vertices.push_back(glm::vec4(kFloorXMax, kFloorY, kFloorZMin, 1.0f));
  floor_vertices.push_back(glm::vec4(kFloorXMin, kFloorY, kFloorZMin, 1.0f));
  floor_faces.push_back(glm::uvec3(0, 1, 2));
  floor_faces.push_back(glm::uvec3(2, 3, 0));

  std::vector<std::string> jpeg_file_names;
  DIR* dir;
  struct dirent* entry;
  CHECK_SUCCESS((dir = opendir("./textures")) != NULL);
  while ((entry = readdir(dir)) != NULL) {
    std::string file_name(entry->d_name);
    std::transform(file_name.begin(), file_name.end(), file_name.begin(),
                   tolower);
    if (file_name.find(".jpg") != std::string::npos) {
      jpeg_file_names.push_back(file_name);
    }
  }
  closedir(dir);

  std::vector<Image> images(jpeg_file_names.size());
  for (int i = 0; i < jpeg_file_names.size(); ++i) {
    std::string file_name = "./textures/" + jpeg_file_names[i];
    LoadJPEG(file_name, &images[i]);
    std::cout << "Loaded '" << file_name << "' width = " << images[i].width
              << " height = " << images[i].height << "\n";
  }



  // Setup the object array object.
  LoadObj("ogre-rigged/ogre.obj", ogre_vertices, ogre_faces);
  LoadBones("ogre-rigged/ogre-skeleton.bf");

  // Setup our VAOs.
  CHECK_GL_ERROR(glGenVertexArrays(kNumVaos, array_objects));

  // Switch to the floor VAO.
  CHECK_GL_ERROR(glBindVertexArray(array_objects[kFloorVao]));

  // Generate buffer objects
  CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kFloorVao][0]));

  // Setup vertex data in a VBO.
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kFloorVao][kVertexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(float) * floor_vertices.size() * 4,
                              &floor_vertices[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(0));

  // Setup element array buffer.
  CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                              buffer_objects[kFloorVao][kIndexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                              sizeof(uint32_t) * floor_faces.size() * 3,
                              &floor_faces[0], GL_STATIC_DRAW));



  // Switch to the skeleton VAO.
  CHECK_GL_ERROR(glBindVertexArray(array_objects[kSkeletonVao]));

  // Generate buffer objects
  CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kSkeletonVao][0]));

  // Setup vertex data in a VBO.
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kSkeletonVao][kVertexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(float) * bone_vertices.size() * 4,
                              &bone_vertices[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(0));

  std::cout << "bone_lines.size: " << bone_lines.size() << std::endl;

  // Setup element array buffer.
  CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                              buffer_objects[kSkeletonVao][kIndexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                              sizeof(uint32_t) * bone_lines.size() * 2,
                             &bone_lines[0], GL_STATIC_DRAW));



    // Switch to the ogre VAO.
  CHECK_GL_ERROR(glBindVertexArray(array_objects[kOgreVao]));

  // Generate buffer objects
  CHECK_GL_ERROR(glGenBuffers(kNumVbos, &buffer_objects[kOgreVao][0]));

  // Setup vertex data in a VBO.
  CHECK_GL_ERROR(
      glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[kOgreVao][kVertexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
                              sizeof(float) * ogre_vertices.size() * 4,
                              &ogre_vertices[0], GL_STATIC_DRAW));
  CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
  CHECK_GL_ERROR(glEnableVertexAttribArray(0));

  // Setup element array buffer.
  CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                              buffer_objects[kOgreVao][kIndexBuffer]));
  CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                              sizeof(uint32_t) * ogre_faces.size() * 3,
                              &ogre_faces[0], GL_STATIC_DRAW));




  // Triangle shaders

  // Setup vertex shader.
  GLuint vertex_shader_id = 0;
  const char* vertex_source_pointer = vertex_shader;
  CHECK_GL_ERROR(vertex_shader_id = glCreateShader(GL_VERTEX_SHADER));
  CHECK_GL_ERROR(
      glShaderSource(vertex_shader_id, 1, &vertex_source_pointer, nullptr));
  glCompileShader(vertex_shader_id);
  CHECK_GL_SHADER_ERROR(vertex_shader_id);

  // Setup geometry shader.
  GLuint geometry_shader_id = 0;
  const char* geometry_source_pointer = geometry_shader;
  CHECK_GL_ERROR(geometry_shader_id = glCreateShader(GL_GEOMETRY_SHADER));
  CHECK_GL_ERROR(
      glShaderSource(geometry_shader_id, 1, &geometry_source_pointer, nullptr));
  glCompileShader(geometry_shader_id);
  CHECK_GL_SHADER_ERROR(geometry_shader_id);

  // Setup floor fragment shader.
  GLuint floor_fragment_shader_id = 0;
  const char* floor_fragment_source_pointer = floor_fragment_shader;
  CHECK_GL_ERROR(floor_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
  CHECK_GL_ERROR(glShaderSource(floor_fragment_shader_id, 1,
                                &floor_fragment_source_pointer, nullptr));
  glCompileShader(floor_fragment_shader_id);
  CHECK_GL_SHADER_ERROR(floor_fragment_shader_id);

  // Let's create our floor program.
  GLuint floor_program_id = 0;
  CHECK_GL_ERROR(floor_program_id = glCreateProgram());
  CHECK_GL_ERROR(glAttachShader(floor_program_id, vertex_shader_id));
  CHECK_GL_ERROR(glAttachShader(floor_program_id, floor_fragment_shader_id));
  CHECK_GL_ERROR(glAttachShader(floor_program_id, geometry_shader_id));

  // Bind attributes.
  CHECK_GL_ERROR(glBindAttribLocation(floor_program_id, 0, "vertex_position"));
  CHECK_GL_ERROR(glBindFragDataLocation(floor_program_id, 0, "fragment_color"));
  glLinkProgram(floor_program_id);
  CHECK_GL_PROGRAM_ERROR(floor_program_id);

  // Get the uniform locations.
  GLint floor_projection_matrix_location = 0;
  CHECK_GL_ERROR(floor_projection_matrix_location =
                     glGetUniformLocation(floor_program_id, "projection"));
  GLint floor_model_matrix_location = 0;
  CHECK_GL_ERROR(floor_model_matrix_location =
                     glGetUniformLocation(floor_program_id, "model"));
  GLint floor_view_matrix_location = 0;
  CHECK_GL_ERROR(floor_view_matrix_location =
                     glGetUniformLocation(floor_program_id, "view"));
  GLint floor_light_position_location = 0;
  CHECK_GL_ERROR(floor_light_position_location =
                     glGetUniformLocation(floor_program_id, "light_position"));

  glm::vec4 light_position = glm::vec4(0.0f, 100.0f, 0.0f, 1.0f);



  // Setup skeleton vertex shader.
  GLuint skel_vertex_shader_id = 0;
  const char* skel_vertex_source_pointer = skeleton_vertex_shader;
  CHECK_GL_ERROR(skel_vertex_shader_id = glCreateShader(GL_VERTEX_SHADER));
  CHECK_GL_ERROR(
      glShaderSource(skel_vertex_shader_id, 1, &skel_vertex_source_pointer, nullptr));
  glCompileShader(skel_vertex_shader_id);
  CHECK_GL_SHADER_ERROR(skel_vertex_shader_id);

  // Setup skeleton geometry shader.
  GLuint skel_geometry_shader_id = 0;
  const char* skel_geometry_source_pointer = skeleton_geometry_shader;
  CHECK_GL_ERROR(skel_geometry_shader_id = glCreateShader(GL_GEOMETRY_SHADER));
  CHECK_GL_ERROR(
      glShaderSource(skel_geometry_shader_id, 1, &skel_geometry_source_pointer, nullptr));
  glCompileShader(skel_geometry_shader_id);
  CHECK_GL_SHADER_ERROR(skel_geometry_shader_id);

  // Setup skeleton fragment shader.
  GLuint skeleton_fragment_shader_id = 0;
  const char* skeleton_fragment_source_pointer = skeleton_fragment_shader;
  CHECK_GL_ERROR(skeleton_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
  CHECK_GL_ERROR(glShaderSource(skeleton_fragment_shader_id, 1,
                                &skeleton_fragment_source_pointer, nullptr));
  glCompileShader(skeleton_fragment_shader_id);
  CHECK_GL_SHADER_ERROR(skeleton_fragment_shader_id);

  // Let's create our skeleton program.
  GLuint skeleton_program_id = 0;
  CHECK_GL_ERROR(skeleton_program_id = glCreateProgram());
  CHECK_GL_ERROR(glAttachShader(skeleton_program_id, skel_vertex_shader_id));
  CHECK_GL_ERROR(glAttachShader(skeleton_program_id, skeleton_fragment_shader_id));
  CHECK_GL_ERROR(glAttachShader(skeleton_program_id, skel_geometry_shader_id));

  // Bind attributes.
  CHECK_GL_ERROR(glBindAttribLocation(skeleton_program_id, 0, "vertex_position"));
  CHECK_GL_ERROR(glBindFragDataLocation(skeleton_program_id, 0, "fragment_color"));
  glLinkProgram(skeleton_program_id);
  CHECK_GL_PROGRAM_ERROR(skeleton_program_id);

  // Get the uniform locations.
  GLint skeleton_projection_matrix_location = 0;
  CHECK_GL_ERROR(skeleton_projection_matrix_location =
                     glGetUniformLocation(skeleton_program_id, "projection"));
  GLint skeleton_model_matrix_location = 0;
  CHECK_GL_ERROR(skeleton_model_matrix_location =
                     glGetUniformLocation(skeleton_program_id, "model"));
  GLint skeleton_view_matrix_location = 0;
  CHECK_GL_ERROR(skeleton_view_matrix_location =
                     glGetUniformLocation(skeleton_program_id, "view"));


  while (!glfwWindowShouldClose(window)) {
    // Setup some basic window stuff.
    glfwGetFramebufferSize(window, &window_width, &window_height);
    glViewport(0, 0, window_width, window_height);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_LESS);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glCullFace(GL_BACK);

    // Compute our view, and projection matrices.
    if (fps_mode)
      center = eye - camera_distance * look;
    else
      eye = center + camera_distance * look;

    view_matrix = glm::lookAt(eye, center, up);
    light_position = glm::vec4(eye, 1.0f);

    aspect = static_cast<float>(window_width) / window_height;
    projection_matrix =
        glm::perspective((float)(kFov * (M_PI / 180.0f)), aspect, kNear, kFar);
    model_matrix = glm::mat4(1.0f);

    // Bind to our floor VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kFloorVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(floor_program_id));

    // Pass uniforms in.
    CHECK_GL_ERROR(glUniformMatrix4fv(floor_projection_matrix_location, 1,
                                      GL_FALSE, &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(floor_model_matrix_location, 1, GL_FALSE,
                                      &floor_model_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(floor_view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));
    CHECK_GL_ERROR(
        glUniform4fv(floor_light_position_location, 1, &light_position[0]));

    // Draw our triangles.
    CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, floor_faces.size() * 3,
                                  GL_UNSIGNED_INT, 0));


    // Bind to our floor VAO.
    CHECK_GL_ERROR(glBindVertexArray(array_objects[kSkeletonVao]));

    // Use our program.
    CHECK_GL_ERROR(glUseProgram(skeleton_program_id));

    // Pass uniforms in.
    CHECK_GL_ERROR(glUniformMatrix4fv(skeleton_projection_matrix_location, 1,
                                      GL_FALSE, &projection_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(skeleton_model_matrix_location, 1, GL_FALSE,
                                      &floor_model_matrix[0][0]));
    CHECK_GL_ERROR(glUniformMatrix4fv(skeleton_view_matrix_location, 1, GL_FALSE,
                                      &view_matrix[0][0]));

    // Draw our triangles.
    CHECK_GL_ERROR(glDrawElements(GL_LINES, bone_lines.size() * 2,
                                  GL_UNSIGNED_INT, 0));



    // // Bind to our floor VAO.
    // CHECK_GL_ERROR(glBindVertexArray(array_objects[kOgreVao]));

    // // Use our program.
    // CHECK_GL_ERROR(glUseProgram(floor_program_id));

    // // Pass uniforms in.
    // CHECK_GL_ERROR(glUniformMatrix4fv(floor_projection_matrix_location, 1,
    //                                   GL_FALSE, &projection_matrix[0][0]));
    // CHECK_GL_ERROR(glUniformMatrix4fv(floor_model_matrix_location, 1, GL_FALSE,
    //                                   &floor_model_matrix[0][0]));
    // CHECK_GL_ERROR(glUniformMatrix4fv(floor_view_matrix_location, 1, GL_FALSE,
    //                                   &view_matrix[0][0]));
    // CHECK_GL_ERROR(
    //     glUniform4fv(floor_light_position_location, 1, &light_position[0]));

    // // Draw our triangles.
    // CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, ogre_faces.size() * 3,
    //                               GL_UNSIGNED_INT, 0));


    // Poll and swap.
    glfwPollEvents();
    glfwSwapBuffers(window);
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  for (int i = 0; i < images.size(); ++i) delete[] images[i].bytes;
  exit(EXIT_SUCCESS);
}
