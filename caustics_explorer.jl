using LinearAlgebra
using StaticArrays
using GLMakie

const TSO3 = SMatrix{3, 3, Float64, 9}
const TCaustic = Tuple{Float64,Float64,TSO3,TSO3,SVector{3,Float64}}
const angres :: Int = 720
const angles :: Vector{Float64} = collect(range(-pi, pi, angres))
const EX = @SVector [1.0, 0.0, 0.0]
const EY = @SVector [0.0, 1.0, 0.0]
const EZ = @SVector [0.0, 0.0, 1.0]

function plane_basis(plane::Symbol)
  if plane == :xy
    return EX, EY
  elseif plane == :xz
    return EX, EZ
  elseif plane == :yz
    return EY, EZ
  else
   error("Unknown plane: $plane")
  end
end

function plane_point(r::Float64, a::Float64,plane::Symbol)
    e1, e2 = plane_basis(plane)
    return r * (cos(a) * e1 + sin(a) * e2)
end



function SO3(yaw::Float64,pitch::Float64,roll::Float64) :: TSO3
  r1 = @SMatrix [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0;0 0 1]
  r2 = @SMatrix [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)]
  r3 = @SMatrix [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)]
  return r1 * r2 * r3
end

function compute_caustic(theta::Float64,phi::Float64,roll::Float64,plane::Symbol,ctype::Symbol) :: TCaustic
  mr = cos(theta)
  pr = sin(phi)
  R1 = SO3(theta,phi,0.0)
  R2 = SO3(theta,phi,roll)
  mp = plane_point(mr,roll,plane)
  #mp = @SVector [mr*cos(roll),0,mr*sin(roll)]
  if ctype == :type1
    return (mr,pr,R1,R2,mp)
  elseif ctype == :type2
    return (pr,mr,R1,R2,mp)
  elseif ctype == :type3
    return (mr,mr,R1,R2,mp)
  elseif ctype == :type4
    return (mr,pr,R2,R1,mp)
  elseif ctype == :type5
    return (pr,mr,R2,R1,mp)
  elseif ctype == :type6
    return (mr,mr,R2,R1,mp)  
  else
   error("Unknown type symbol")
  end
end

function caustic_pt(c::TCaustic,ca::Float64,plane::Symbol) :: SVector{3,Float64}
  mr = c[1]
  pr = c[2]
  R1 = c[3]
  R2 = c[4]
  mp = c[5]
  p = plane_point(pr,ca,plane) 
  p1 = R1 * p
  p2 = R2 * p1
  return mp + p2
end

function compute_pts(theta::Float64, phi::Float64, mid_plane::Symbol,point_plane::Symbol,ctype::Symbol)
  pts::Vector{SVector{3,Float64}} = []
  phi_ = ctype in [:type3,:type6] ? 0.0 : phi 
  for roll in angles
    c = compute_caustic(theta,phi_,roll,mid_plane,ctype)
    for a in angles
      push!(pts,caustic_pt(c,a,point_plane))
    end
  end
  return pts
end

fig = Figure(size = (1600, 850))
ax=Axis3(fig[1:10, 1:10])
theta_slider = Slider(fig[11,1:3],range = range(-pi, pi,length=1000),startvalue=0.7)
Label(fig[12,1:3],lift(theta -> "theta:$theta",theta_slider.value),fontsize = 18)    
phi_slider = Slider(fig[11, 4:6],range = range(0, pi, length=1000),startvalue=1.2)
Label(fig[12, 4:6],lift(phi -> "phi:$phi",phi_slider.value),fontsize=18)
mp_menu = Menu(fig[11, 7:8],options = [("XY", :xy),("XZ", :xz),("YZ", :yz)],default = "XZ",width = 100)
Label(fig[12, 7:8], "mid plane")
pp_menu = Menu(fig[11, 9:10],options = [("XY", :xy),("XZ", :xz),("YZ", :yz)],default ="XY",width = 100)
Label(fig[12, 9:10], "point plane")
ct_menu = Menu(fig[11, 11:12],options = [("type1", :type1),("type2", :type2),("type3", :type3),("type4", :type4),("type5", :type5),("type6", :type6)],default ="type1",width = 100)
Label(fig[12, 11:12], "caustic type")


pts = Observable(compute_pts(0.7,1.2,mp_menu.selection[],pp_menu.selection[],ct_menu.selection[]))
lines!(ax, pts, linewidth=0.5)


onany(theta_slider.value, phi_slider.value,mp_menu.selection,pp_menu.selection,ct_menu.selection) do theta,phi,mp,pp,ct
  #print(ct == :type1)
  #print(ct == :type2)
  pts[] = compute_pts(theta,phi,mp,pp,ct)
end

display(fig)
