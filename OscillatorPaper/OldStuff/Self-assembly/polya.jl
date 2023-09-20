using Plots

function random_walk_2d_2points_trajectory(steps)
    # create a 2D matrix of size 50x50
    lattice = zeros(10, 10)
    # starting positions
    x1, y1 = rand(1:10), rand(1:10)
    x2, y2 = rand(1:10), rand(1:10)
    lattice[x1, y1] = 1
    lattice[x2, y2] = 2
    # create a matrix to store the trajectories
    trajectory1 = [(x1, y1)]
    trajectory2 = [(x2, y2)]
    
    for i in 1:steps
        # randomly choose new position for point 1
        while true
            dx1, dy1 = rand([-1, 0, 1]), rand([-1, 0, 1])
            x1_new, y1_new = x1 + dx1, y1 + dy1
            if x1_new >= 1 && x1_new <= size(lattice, 1) && y1_new >= 1 && y1_new <= size(lattice, 2)
                x1, y1 = x1_new, y1_new
                break
            end
        end
        
        # randomly choose new position for point 2
        while true
            dx2, dy2 = rand([-1, 0, 1]), rand([-1, 0, 1])
            x2_new, y2_new = x2 + dx2, y2 + dy2
            if x2_new >= 1 && x2_new <= size(lattice, 1) && y2_new >= 1 && y2_new <= size(lattice, 2)
                x2, y2 = x2_new, y2_new
                break
            end
        end
        
        # check if the two point are in the same position
        if x1 == x2 && y1 == y2
            return lattice, trajectory1, trajectory2
        end
        
        lattice[x1, y1] = 1
        lattice[x2, y2] = 2
        push!(trajectory1, (x1, y1))
        push!(trajectory2, (x2, y2))
    end
    return lattice, trajectory1, trajectory2
end



random_walk_2d_2points_trajectory(100)


function random_walk_3d_2points_trajectory(steps)
    # create a 3D matrix
    lattice = zeros(10, 10, 10)
    # starting positions
    x1, y1, z1 = rand(1:10), rand(1:10), rand(1:10)
    x2, y2, z2 = rand(1:10), rand(1:10), rand(1:10)
    lattice[x1, y1, z1] = 1
    lattice[x2, y2, z2] = 2
    # create a matrix to store the trajectories
    trajectory1 = [(x1, y1, z1)]
    trajectory2 = [(x2, y2, z2)]
    for i in 1:steps
        # randomly choose new position for point 1
        while true
            dx1, dy1, dz1 = rand([-1, 0, 1]), rand([-1, 0, 1]), rand([-1, 0, 1])
            x1_new, y1_new, z1_new = x1 + dx1, y1 + dy1, z1 + dz1
            if x1_new >= 1 && x1_new <= size(lattice, 1) && y1_new >= 1 && y1_new <= size(lattice, 2) && z1_new >= 1 && z1_new <= size(lattice, 3)
                x1, y1, z1 = x1_new, y1_new, z1_new
                break
            end
        end

        while true
            dx2, dy2, dz2 = rand([-1, 0, 1]), rand([-1, 0, 1]), rand([-1, 0, 1])
            x2_new, y2_new, z2_new = x2 + dx2, y2 + dy2, z2 + dz2
            if x2_new >= 1 && x2_new <= size(lattice, 1) && y2_new >= 1 && y2_new <= size(lattice, 2) && z2_new >= 1 && z2_new <= size(lattice, 3)
                x2, y2, z2 = x2_new, y2_new, z2_new
                break
            end
        end

        # check if the two points are in the same position
        if x1 == x2 && y1 == y2 && z1 == z2
            return lattice, trajectory1, trajectory2
        end
        
        lattice[x1, y1, z1] = 1
        lattice[x2, y2, z2] = 2
        push!(trajectory1, (x1, y1, z1))
        push!(trajectory2, (x2, y2, z2))
    end
    return lattice, trajectory1, trajectory2
end

random_walk_3d_2points_trajectory(100)




default(dpi = 200, linewidth = 2, markersize = 5, legend = :topright, size = (500, 500))

# number of steps for the random walk
steps = 500

# generate the 2D random walk trajectory
_, trajectory1, trajectory2 = random_walk_2d_2points_trajectory(steps)

# check if the two points co-localize
if size(trajectory1, 1) < steps
    scatter!([trajectory1[end]], color = :red)
end

# create an animation for the 2D random walk trajectory
animation1 = @animate for i in 1:steps
    if i == size(trajectory1, 1)
        scatter!([trajectory1[end]], color = :green, markersize = 10, label = "Reaction!")
    elseif i < size(trajectory1, 1)
        # @info "Plotting the $i th step"
        plot(trajectory1[1:i], title = "2D Search", label = "Point 1", arrow = true, xlims = (1,10), ylims = (1,10))
        plot!(trajectory2[1:i], label = "Point 2", arrow = true)
    end
end

gif(animation1, "2d_random_walk.gif")

# repeat the same for 3D random walk trajectory
_, trajectory3, trajectory4 = random_walk_3d_2points_trajectory(steps)
animation2 = @animate for i in 1:steps
    if i == size(trajectory3, 1)
        scatter!([trajectory3[end]], color = :green, markersize = 10, label = "Reaction!")
    elseif i < size(trajectory3, 1)
        # @info "Plotting the $i th step"
        plot(trajectory3[1:i], title = "3D Search", label = "Point 1", arrow = true, xlims = (1,10), ylims = (1,10), zlims = (1,10))
        plot!(trajectory4[1:i], label = "Point 2", arrow = true)
    end
end

gif(animation2, "3d_random_walk.gif")

