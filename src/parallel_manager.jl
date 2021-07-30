# 管理不同worker上的任务，并给出进度
function auto_manger(task)
    # 总共的任务数量
    num_of_task = length(task)
    # println(num_of_task)
    # println(task[num_of_task])
    # 获取远程进程数量
    numof_workers = length(workers())
    # println(numof_workers)
    result = give_empty_dataframe_for_format()
    workersFuture = [Future(1) for i in 1:numof_workers]
    # 记录完成任务ID在列表中的位置
    c = Channel{Int64}(numof_workers)

    progress = 0.0

    # 同时启动numof_workers数量的进程
    for i in 1:numof_workers
        # println("$num_of_task")
        # @async put!(workersFuture[i], remotecall_fetch(cal_one_task, workers()[i], task[i])) & put!(c, i)
        @async begin; put!(workersFuture[i], remotecall_fetch(cal_one_task, workers()[i], task[i])); put!(c, i); end
        # println("first")
    end
    # print("？")
    # i记录了下次数据在result中存放的位置
    # (i - 1)即为已经记录的数据数量
    # 并且(i - 1 + numof_workers)为已经启动的进程数
    # 当已经启动的进程数量小于num_of_task时，不断记录已经完成任务的进程的结果，并在完成任务的进程再启动完成任务
    # 当启动num_of_task数量的进程后，等待剩余正在运行的进程并返回结果
    i = 1
    # while (i - 1 + numof_workers) < num_of_task
    test_num = 0
    while  true
            
        # yield()
        # 完成任务的进程序号
        j = take!(c)
        # 取出结果
        append!(result, fetch(workersFuture[j]))
        # println(result[i])
        workersFuture[j] = Future(1)
        # println("up $i #")
        # println(task[i + numof_workers].W)
        # println(workers()[j])
        # println((i - 1 + numof_workers) < num_of_task)
        # print(i + numof_workers)
        # print(", j $j \n")

        # 未启动过num_of_task的进程则继续启动
        if (i - 1 + numof_workers) < num_of_task
            # @async put!(workersFuture[j], remotecall_fetch(cal_one_task, workers()[j], task[i + numof_workers])) & put!(c, j)
            the_task = task[i + numof_workers]
            the_id::Int = workers()[j]
            the_j = j
            # k = i + numof_workers
            # println(the_task)
            # println(the_id)
            # @async begin; println(i); put!(workersFuture[j], remotecall_fetch(cal_one_task, workers()[j], task[i + numof_workers])); put!(c, j); end
            @async begin; put!(workersFuture[the_j], remotecall_fetch(cal_one_task, the_id, the_task)); put!(c, the_j); end
            test_num = test_num + 1
            # println(i)
            # println("test_num $test_num")
            # 进度条
            if (i - 1)/num_of_task > 0.001 + progress
                progress = i/num_of_task
                print("\r")
                print(round(progress*100, digits = 1))
                print("%")
            end
            # 进度条
        # 未记录num_of_task数量的数量则继续记录
        elseif i < num_of_task
            nothing
        else
            break
        end
        # println(length(result))
        # yield()
        i = i + 1



    end

    close(c)
    print("\r100.0%")
    print("  all done @")
    println(now())
    result

end
