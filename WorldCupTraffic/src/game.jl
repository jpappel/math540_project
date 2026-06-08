using Distributions

function SoccerGame(t1)
    time = 0;
    score = [0 0];

    # probabilities of each team scoring
    t2 = 1-t1;

    # the main game
    while time < 90
        time += rand(Exponential(lambda));
        if time < 90
            g = rand(Uniform(0,1))
            if g < t1
                score[1] += 1
            else
                score[2] += 1
            end
        end
    end
    #overtime
    if score[1] == score[2]
        while time < 120
        time += rand(Exponential(lambda));
            if time < 120
                g = rand(Uniform(0,1))
                if g < t1
                    score[1] += 1
                else
                    score[2] += 1
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
            end
            if gt2 < t2
                score[2] += 1
            end
        end
    end
    return score
end

export SoccerGame
