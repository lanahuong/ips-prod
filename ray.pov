camera
{
  location <5,5,5>
  look_at <0,0,0>
  angle 0
}


box
{
 <0,0,0><1,1,1>
 hollow
 pigment
 { rgbf 1 }
 interior
 {
  media
  {
   emission 0.7
   absorption 1.0
   intervals 10
   density
   {
    density_file df3 "example.df3"
    scale<0.7, 0.7, 0.7>
    interpolate 1
    color_map
    {
     [0.0   rgb <0.0,0.0,0.0>]
     [0.1   rgb <0.0,0.0,2.0>]
     [1.0   rgb <1.0,1.0,1.0>]
    }
   }
  }
 }
 scale <10,10,20>
 translate <-5,-5,-10>
 rotate y*90
}