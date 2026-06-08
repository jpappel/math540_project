using Distributions

function SoccerGame(t1)
    goals = []
    time = 0;
    score = [0 0];
    lambda = 36;

    # probabilities of each team scoring
    t2 = 1-t1;

    # the main game
    while time < 90
        time += rand(Exponential(lambda));
        if time < 90
            g = rand(Uniform(0,1))
            if g < t1
                score[1] += 1
                push!(goals,(time,1))
            else
                score[2] += 1
                push!(goals,(time,2))
            end
        end
    end
    #overtime
    if score[1] == score[2]
        time = 90;
        while time < 120
            time += rand(Exponential(lambda));
            if time < 120
                g = rand(Uniform(0,1))
                if g < t1
                    score[1] += 1
                push!(goals,(time,1))
                else
                    score[2] += 1
                    push!(goals,(time,2))
                end
            end
        end
    end
    #penalty shootout
    if score[1] == score[2]
        while score[1] == score[2]
            gt1 = rand(Uniform(0,1))
            gt2 = rand(Uniform(0,1))
            if gt1 < t1
                score[1] += 1
                push!(goals,("shootout",1))
            end
            if gt2 < t2
                score[2] += 1
                push!(goals,("shootout",2))
            end
        end
    end
    return score, goals
end

export SoccerGame
